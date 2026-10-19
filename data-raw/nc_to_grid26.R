#!/usr/bin/env Rscript

# Convert a monthly downscaled-climate NetCDF (e.g. thetao or o2, one depth
# layer, one scenario) to an sf object on the pacea grid26 (2 km inshore,
# 6 km offshore, BC Albers), in the same format as the pacea BCCM outputs:
# a `pacea_st` sf tibble with one column per month named "YYYY_M" and a
# "units" attribute.
#
# Each grid26 cell gets the area-weighted mean of the NetCDF cells it overlaps,
# ignoring NA (land) cells, so coastal cells partly over land still get a value.
# grid26 cells with no data (e.g. outside the NetCDF domain) are dropped.
#
# Can be run from the command line (see usage()) or source()d and called via
# nc_to_grid26().

suppressPackageStartupMessages({
  library(ncdf4)
  library(terra)
  library(sf)
  library(Matrix)
})

parse_args <- function(args) {
  defaults <- list(
    input = NULL,
    output = NULL,
    var = NULL,
    percentile = "50",
    years = "all",
    chunk_months = "60",
    digits = "6",
    overwrite = FALSE
  )

  for (arg in args) {
    if (!startsWith(arg, "--")) {
      stop("All arguments must use the --name=value form.")
    }

    if (identical(arg, "--overwrite")) {
      defaults$overwrite <- TRUE
      next
    }

    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(pieces) != 2L) {
      stop("Argument must be in --name=value form: ", arg)
    }

    key <- gsub("-", "_", pieces[[1]])
    if (!key %in% names(defaults)) {
      stop("Unknown argument: --", pieces[[1]])
    }

    defaults[[key]] <- pieces[[2]]
  }

  if (is.null(defaults$input) || !nzchar(defaults$input)) {
    stop("Missing required argument: --input=/path/to/file.nc")
  }

  defaults$percentile <- as.numeric(defaults$percentile)
  defaults$chunk_months <- as.integer(defaults$chunk_months)
  defaults$digits <- as.integer(defaults$digits)
  defaults$years <- parse_years(defaults$years)
  defaults$overwrite <- isTRUE(defaults$overwrite) ||
    identical(tolower(defaults$overwrite), "true")
  defaults
}

usage <- function() {
  cat(
    paste(
      "Usage:",
      "Rscript data-raw/nc_to_grid26.R \\",
      "  --input=path/to/file.nc \\",
      "  [--output=path/to/output.rds] \\",
      "  [--var=thetao]  (default: the only time/lat/lon variable, or *_mean) \\",
      "  [--percentile=50]  (ignored if the file has no percentiles dimension) \\",
      "  [--years=all|1993:2024,2030,2040] \\",
      "  [--chunk-months=60] \\",
      "  [--digits=6]  (decimal places kept; fewer = much smaller files) \\",
      "  [--overwrite]",
      sep = "\n"
    ),
    "\n"
  )
}

# "1993:2024,2030,2040" -> c(1993:2024, 2030, 2040); "all" -> NULL
parse_years <- function(x) {
  if (is.null(x) || identical(x, "all")) {
    return(NULL)
  }
  pieces <- strsplit(x, ",", fixed = TRUE)[[1]]
  out <- unlist(lapply(pieces, function(p) {
    ends <- as.integer(strsplit(p, ":", fixed = TRUE)[[1]])
    if (anyNA(ends) || !length(ends) %in% 1:2) {
      stop("Could not parse --years element: ", p)
    }
    seq(ends[[1]], ends[[length(ends)]])
  }))
  sort(unique(out))
}

detect_main_var <- function(nc) {
  candidates <- names(nc$var)[vapply(
    nc$var,
    function(x) {
      dims <- vapply(x$dim, `[[`, character(1), "name")
      all(c("time", "lat", "lon") %in% dims)
    },
    logical(1)
  )]

  # "stats" files hold <var>_mean, _min, _max, _stdev; default to the mean
  if (length(candidates) > 1L) {
    candidates <- grep("_mean$", candidates, value = TRUE)
  }
  if (length(candidates) != 1L) {
    stop("Could not pick a variable with time/lat/lon dimensions. Use --var.")
  }
  candidates
}

parse_nc_time <- function(time_values, units_string) {
  pieces <- strsplit(units_string, " since ", fixed = TRUE)[[1]]
  if (length(pieces) != 2L) {
    stop("Unsupported time units: ", units_string)
  }

  origin <- as.POSIXct(trimws(pieces[[2]]), tz = "UTC")
  if (is.na(origin)) {
    stop("Could not parse NetCDF time origin: ", pieces[[2]])
  }

  multiplier <- switch(
    trimws(pieces[[1]]),
    days = , day = 86400,
    hours = , hour = 3600,
    minutes = , minute = 60,
    seconds = , second = 1,
    stop("Unsupported time unit: ", pieces[[1]])
  )

  as.Date(origin + time_values * multiplier, tz = "UTC")
}

# Label stored in attr(x, "units"), matching pacea's style where we know it
units_label <- function(nc, var_name) {
  long_name <- ncatt_get(nc, var_name, "long_name")
  units <- ncatt_get(nc, var_name, "units")
  if (units$hasatt && units$value %in% c("degC", "degrees_C", "Celsius")) {
    return("Temperature (°C)")
  }
  paste0(
    if (long_name$hasatt) long_name$value else var_name,
    if (units$hasatt) paste0(" (", units$value, ")") else ""
  )
}

# Check the lon/lat axes are regularly spaced and return their cell edges
regular_axis_edges <- function(x, name) {
  d <- diff(x)
  if (any(abs(d - mean(d)) > 1e-3 * abs(mean(d)))) {
    stop(name, " is not regularly spaced; this script assumes a regular grid.")
  }
  step <- abs(mean(d))
  c(min(x) - step / 2, max(x) + step / 2)
}

# Sparse matrix W (n_polygons x n_nc_cells) of overlap fractions between each
# target polygon and each NetCDF cell. NetCDF cells are indexed in the order
# ncvar_get() returns them for (lon, lat): lon varies fastest.
build_weights <- function(target, lon, lat) {
  if (!all(diff(lon) > 0)) {
    stop("Longitude must be increasing.")
  }
  lat_increasing <- all(diff(lat) > 0)
  if (!lat_increasing && !all(diff(lat) < 0)) {
    stop("Latitude must be monotonic.")
  }

  nlon <- length(lon)
  nlat <- length(lat)
  x_edges <- regular_axis_edges(lon, "Longitude")
  y_edges <- regular_axis_edges(lat, "Latitude")

  template <- rast(
    nrows = nlat,
    ncols = nlon,
    extent = ext(x_edges[[1]], x_edges[[2]], y_edges[[1]], y_edges[[2]]),
    crs = "EPSG:4326"
  )
  values(template) <- seq_len(ncell(template))

  target_ll <- vect(st_transform(st_geometry(target), "EPSG:4326"))
  e <- terra::extract(template, target_ll, cells = TRUE, exact = TRUE, ID = TRUE)
  e <- e[!is.na(e$cell) & e$fraction > 0, ]

  # terra cells run row-wise from the north; convert to NetCDF (lon, lat) order
  row <- (e$cell - 1L) %/% nlon + 1L
  col <- (e$cell - 1L) %% nlon + 1L
  lat_idx <- if (lat_increasing) nlat - row + 1L else row
  nc_cell <- col + (lat_idx - 1L) * nlon

  sparseMatrix(
    i = e$ID,
    j = nc_cell,
    x = e$fraction,
    dims = c(nrow(target), nlon * nlat)
  )
}

# Read a contiguous run of time steps as a (n_cells x n_time) matrix
read_time_run <- function(nc, var_name, dim_names, time_idx, percentile_idx) {
  start <- rep(1L, length(dim_names))
  count <- rep(1L, length(dim_names))
  for (i in seq_along(dim_names)) {
    switch(
      dim_names[[i]],
      lon = , lat = {
        count[[i]] <- -1L
      },
      time = {
        start[[i]] <- time_idx[[1]]
        count[[i]] <- length(time_idx)
      },
      percentiles = {
        start[[i]] <- percentile_idx
      }
    )
  }

  raw <- ncvar_get(nc, var_name, start = start, count = count, collapse_degen = FALSE)
  keep <- match(c("lon", "lat", "time"), dim_names)
  raw <- aperm(raw, c(keep, setdiff(seq_along(dim_names), keep)))
  matrix(raw, ncol = length(time_idx))
}

# Area-weighted mean of the non-NA cells under each polygon
weighted_means <- function(W, vals) {
  ok <- !is.na(vals)
  vals[!ok] <- 0
  num <- as.matrix(W %*% vals)
  den <- as.matrix(W %*% ok)
  out <- num / den
  out[den == 0] <- NA_real_
  out
}

nc_to_grid26 <- function(
  input,
  output = NULL,
  var = NULL,
  percentile = 50,
  years = NULL,
  chunk_months = 60L,
  digits = 6L,
  overwrite = FALSE,
  grid = pacea::grid26
) {
  if (!file.exists(input)) {
    stop("Input file does not exist: ", input)
  }

  nc <- nc_open(input)
  on.exit(nc_close(nc), add = TRUE)

  var_name <- if (is.null(var)) detect_main_var(nc) else var
  dim_names <- vapply(nc$var[[var_name]]$dim, `[[`, character(1), "name")

  percentile_idx <- NULL
  if ("percentiles" %in% dim_names) {
    percentile_values <- ncvar_get(nc, "percentiles")
    percentile_idx <- which(abs(percentile_values - percentile) < 1e-8)
    if (length(percentile_idx) != 1L) {
      stop(
        "Percentile ", percentile, " not found. Available: ",
        paste(percentile_values, collapse = ", ")
      )
    }
  }

  if (is.null(output)) {
    stem <- tools::file_path_sans_ext(basename(input))
    # e.g. "_p50" for percentile files, "_mean" for thetao_mean in stats files
    suffix <- if (is.null(percentile_idx)) {
      sub("^[^_]+", "", var_name)
    } else {
      paste0("_p", gsub("\\.", "p", percentile))
    }
    output <- file.path(dirname(input), paste0(stem, suffix, "_grid26.rds"))
  }
  if (file.exists(output) && !isTRUE(overwrite)) {
    stop("Output exists already. Use overwrite = TRUE to replace: ", output)
  }

  dates <- parse_nc_time(ncvar_get(nc, "time"), ncatt_get(nc, "time", "units")$value)
  year <- as.integer(format(dates, "%Y"))
  month <- as.integer(format(dates, "%m"))
  if (anyDuplicated(paste(year, month))) {
    stop("Time axis has more than one step per month; expected monthly data.")
  }

  time_idx <- if (is.null(years)) seq_along(dates) else which(year %in% years)
  if (length(time_idx) == 0L) {
    stop("No time steps selected.")
  }

  lon <- as.vector(ncvar_get(nc, "lon"))
  lat <- as.vector(ncvar_get(nc, "lat"))

  message("Input: ", normalizePath(input))
  message("Variable: ", var_name, if (!is.null(percentile_idx)) paste0(" (p", percentile, ")"))
  message(
    "Months: ", length(time_idx), " of ", length(dates),
    " (", format(dates[min(time_idx)], "%Y-%m"), " to ", format(dates[max(time_idx)], "%Y-%m"), ")"
  )

  W <- build_weights(grid, lon, lat)

  # Read in contiguous runs of at most chunk_months time steps
  run_id <- cumsum(c(1L, diff(time_idx) != 1L))
  chunks <- unlist(
    lapply(split(time_idx, run_id), function(x) split(x, ceiling(seq_along(x) / chunk_months))),
    recursive = FALSE
  )

  out <- matrix(NA_real_, nrow = nrow(grid), ncol = length(time_idx))
  col_offset <- 0L
  for (chunk in chunks) {
    vals <- read_time_run(nc, var_name, dim_names, chunk, percentile_idx)
    cols <- col_offset + seq_along(chunk)
    out[, cols] <- weighted_means(W, vals)
    col_offset <- col_offset + length(chunk)
    message("  ", format(dates[max(chunk)], "%Y-%m"), " done")
  }

  colnames(out) <- paste(year[time_idx], month[time_idx], sep = "_")
  has_data <- rowSums(!is.na(out)) > 0L
  message("grid26 cells with data: ", sum(has_data), " of ", nrow(grid))

  result <- sf::st_as_sf(
    tibble::as_tibble(round(out[has_data, , drop = FALSE], digits = digits)),
    geometry = st_geometry(grid)[has_data]
  )
  class(result) <- c("pacea_st", "sf", "tbl_df", "tbl", "data.frame")
  attr(result, "units") <- units_label(nc, var_name)

  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, output, compress = "xz")
  message("Wrote: ", normalizePath(output))
  invisible(output)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0L || any(args %in% c("--help", "-h"))) {
    usage()
    return(invisible(NULL))
  }

  opts <- parse_args(args)
  nc_to_grid26(
    input = opts$input,
    output = opts$output,
    var = opts$var,
    percentile = opts$percentile,
    years = opts$years,
    chunk_months = opts$chunk_months,
    digits = opts$digits,
    overwrite = opts$overwrite
  )
}

if (sys.nframe() == 0L) {
  main()
}
