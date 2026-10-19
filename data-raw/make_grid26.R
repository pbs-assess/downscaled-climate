# Convert every downscaled NetCDF on the external drive to a pacea grid26 .rds.
#
# Expected layout (unzip ssp245.zip and ssp585.zip next to ssp126/ first):
#   <raw_dir>/<scenario>/ensemble/<var>/<var>_ensemble_sd+ba_<layer>_{perc,stats}_<scenario>.nc
# e.g. ssp126/ensemble/thetao/thetao_ensemble_sd+ba_bottom_stats_ssp126.nc
#      ssp126/ensemble/tos/tos_ensemble_sd+ba_surface_depth_0.5_perc_ssp126.nc
#
# Run from the project root: Rscript data-raw/make_grid26.R

source("data-raw/nc_to_grid26.R")

raw_dir <- "/Volumes/Toshiba-Storage/downscaling"
output_dir <- "output/grid26"

# "mean" uses the *_stats_* files (ensemble mean); "median" the *_perc_* files (p50)
statistic <- "mean"

# All historical years, then one year per future decade (data end in 2099).
# NULL keeps every year.
years <- c(1993:2025, seq(2030, 2090, by = 10))

# Output names follow pacea's bccm_<layer>_<variable> style, e.g.
# thetao_ensemble_sd+ba_bottom_stats_ssp126.nc -> bottom_temperature_ssp126_mean.rds
variable_names <- c(thetao = "temperature", tos = "temperature")

# Decimal places kept per variable. File size scales roughly with this:
# for temperature, 6 dp -> ~150 MB, 2 dp (0.01 degC) -> ~40 MB per file.
# Check the units before choosing (e.g. oxygen in mol m-3 needs more digits).
variable_digits <- c(thetao = 2L, tos = 2L)

output_name <- function(nc_file, statistic) {
  stem <- tools::file_path_sans_ext(basename(nc_file))
  var <- sub("_.*", "", stem)
  if (!var %in% names(variable_names)) {
    stop("Add '", var, "' to variable_names in make_grid26.R")
  }
  layer <- regmatches(stem, regexpr("(bottom|surface)", stem))
  scenario <- regmatches(stem, regexpr("ssp[0-9]+", stem))
  if (length(layer) != 1L || length(scenario) != 1L) {
    stop("Could not parse layer/scenario from file name: ", basename(nc_file))
  }
  paste0(layer, "_", variable_names[[var]], "_", scenario, "_", statistic, ".rds")
}

file_type <- switch(statistic, mean = "stats", median = "perc")
files <- list.files(
  raw_dir,
  pattern = paste0("_", file_type, "_ssp[0-9]+\\.nc$"),
  recursive = TRUE,
  full.names = TRUE
)
if (length(files) == 0L) {
  stop("No *_", file_type, "_*.nc files found under ", raw_dir)
}

for (f in files) {
  out <- file.path(output_dir, output_name(f, statistic))
  if (file.exists(out)) {
    message("Skipping existing: ", out)
    next
  }
  var <- sub("_.*", "", basename(f))
  nc_to_grid26(
    f,
    output = out,
    percentile = 50,
    years = years,
    digits = variable_digits[[var]]
  )
}
