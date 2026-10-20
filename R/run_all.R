# R/run_all.R -- execute the full pipeline from the project root
# Usage: Rscript R/run_all.R
# Scripts 05 and 06 take several minutes.

scripts <- c(
  "R/01_prepare_data.R",
  "R/02_eda.R",
  "R/03_unmet_demand.R",
  "R/04_model.R",
  "R/05_monte_carlo.R",
  "R/05b_block_maxima.R",
  "R/06_calibration.R"
)

for (s in scripts) {
  cat("\n################ ", s, " ################\n")
  t0 <- Sys.time()
  source(s, echo = FALSE, keep.source = TRUE)
  cat("Elapsed:", format(Sys.time() - t0), "\n")
}

cat("\nDone. Figures in figures/, tables in outputs/.\n")
