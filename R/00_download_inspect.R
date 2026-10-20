# R/00_download_inspect.R
# Run this FIRST. It downloads the Zenodo files we need and prints their structure.
# Paste the console output back so column names/patterns in 01_prepare_data.R can be fixed.

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
options(timeout = 600)

base  <- "https://zenodo.org/records/14983362/files/"
files <- c("corrected_daily_energy_met_MU.csv",
           "corrected_peak_met_MW.csv",
           "peak_unmet_MW.csv",
           "Tamil_Nadu.csv")

for (f in files) {
  dest <- file.path("data/raw", f)
  if (!file.exists(dest)) {
    download.file(paste0(base, f, "?download=1"), dest, mode = "wb")
  }
}

for (f in files) {
  cat("\n=====", f, "=====\n")
  x <- read.csv(file.path("data/raw", f), check.names = FALSE)
  cat("dim:", dim(x), "\n")
  cat("first column:", names(x)[1], "| first/last value:",
      as.character(x[1, 1]), "/", as.character(x[nrow(x), 1]), "\n")
  cat("columns (first 60):\n"); print(head(names(x), 60))
  cat("head:\n"); print(head(x[, seq_len(min(6, ncol(x)))], 3))
  tn <- grep("tamil", names(x), ignore.case = TRUE, value = TRUE)
  if (length(tn)) {
    cat("Tamil Nadu column:", tn, "| NAs:", sum(is.na(x[[tn[1]]])),
        "| summary:\n"); print(summary(x[[tn[1]]]))
  }
}
