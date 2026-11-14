# Requires completed R_results/analysis_objects.rds and native model tables.
# This also fits the additional whole-metaproteome exclusion model for Supp-F8.
local({
  previous <- getOption("lobaria.run_mode")
  on.exit(options(lobaria.run_mode = previous), add = TRUE)
  options(lobaria.run_mode = "figures")
  source("Lobaria_complete_JPG_only.R", local = globalenv())
})
