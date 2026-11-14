# Full analysis. Run from the repository root.
local({
  previous <- getOption("lobaria.run_mode")
  on.exit(options(lobaria.run_mode = previous), add = TRUE)
  options(lobaria.run_mode = "full")
  source("Lobaria_complete_JPG_only.R", local = globalenv())
})
