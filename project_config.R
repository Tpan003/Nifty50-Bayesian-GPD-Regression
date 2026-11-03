# Project-wide settings. Keep modelling choices in one place so that the
# analysis, diagnostics, and report use the same values.

project_config <- function() {
  list(
    threshold = 2,                       # Daily loss threshold, in percentage points
    rolling_observations = 1250L,       # Number of prior Nifty trading dates
    minimum_training_exceedances = 40L,
    forecast_step = 1L,
    seed = 20261018L,
    prior_specification = gpd_prior_specification()
  )
}
