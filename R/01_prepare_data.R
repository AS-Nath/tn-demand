# R/01_prepare_data.R -- build one clean daily table for Tamil Nadu
source("R/specs.R")

raw <- "data/raw"
UNMET_NA_AS_ZERO <- FALSE   # TRUE only if you decide missing unmet == no deficit

# Regex patterns for weather columns in Tamil_Nadu.csv.
# Exact column names confirmed from 00_download_inspect.R output.
# Each pattern must match exactly one column.
WX_PATTERNS <- c(
  tmax  = "^2m_temperature_max$",
  tmean = "^2m_temperature_mean$",
  dew   = "^2m_dewpoint_temperature_mean$",
  solar = "^surface_solar_radiation_downwards_mean$",
  cloud = "^total_cloud_cover_mean$"
)
WX_OPTIONAL <- c(utci = "^utci_mean$")   # used in EDA only

# ---- helpers -----------------------------------------------------------------
parse_dates <- function(v) {
  if (inherits(v, "Date")) return(v)
  if (inherits(v, "POSIXct")) return(as.Date(v))
  as.Date(lubridate::parse_date_time(as.character(v),
                                     orders = c("Ymd HMS", "Ymd", "dmY HMS", "dmY", "mdY"), quiet = TRUE))
}

read_dated <- function(f) {
  x <- readr::read_csv(file.path(raw, f), show_col_types = FALSE)
  names(x)[1] <- "date"
  x$date <- parse_dates(x$date)
  if (anyNA(x$date)) stop("Unparseable dates in ", f)
  x
}

pick <- function(x, pattern, required = TRUE) {
  nm <- grep(pattern, names(x), ignore.case = TRUE, value = TRUE)
  if (length(nm) == 1) return(nm)
  msg <- sprintf("Pattern '%s' matched %d columns: [%s]\nAll columns: %s",
                 pattern, length(nm), paste(nm, collapse = ", "),
                 paste(names(x), collapse = ", "))
  if (required) stop(msg, call. = FALSE)
  message(msg); NA_character_
}

tn_series <- function(f, newname) {
  x   <- read_dated(f)
  nm  <- pick(x, "tamil")
  out <- tibble(date = x$date, v = as.numeric(x[[nm]]))
  
  # One row per date. Duplicates that agree collapse to that value; duplicates that
  # conflict are set to NA (no way to know which scrape is right).
  out <- out %>%
    group_by(date) %>%
    summarise(conflict = n_distinct(v, na.rm = TRUE) > 1,
              v = if (n_distinct(v, na.rm = TRUE) > 1) NA_real_ else v[!is.na(v)][1],
              .groups = "drop")
  if (any(out$conflict))
    message(f, ": conflicting duplicates set to NA on ",
            paste(format(out$date[out$conflict]), collapse = ", "))
  out %>% select(-conflict) %>% rename(!!newname := v)
}

# ---- demand side ---------------------------------------------------------------
energy    <- tn_series("corrected_daily_energy_met_MU.csv", "energy_MU")
peak_met  <- tn_series("corrected_peak_met_MW.csv", "peak_met")
peak_unm  <- tn_series("peak_unmet_MW.csv", "peak_unmet")

demand <- reduce(list(energy, peak_met, peak_unm), full_join, by = "date")

# ---- weather side --------------------------------------------------------------
wx_raw <- read_dated("Tamil_Nadu.csv")
if (anyDuplicated(wx_raw$date))
  stop("Tamil_Nadu.csv has several rows per date (sub-daily?). Aggregate first.")

wx <- tibble(date = wx_raw$date)
for (v in names(WX_PATTERNS)) wx[[v]] <- as.numeric(wx_raw[[pick(wx_raw, WX_PATTERNS[[v]])]])
for (v in names(WX_OPTIONAL)) {
  nm <- pick(wx_raw, WX_OPTIONAL[[v]], required = FALSE)
  if (!is.na(nm)) wx[[v]] <- as.numeric(wx_raw[[nm]])
}

# ERA5-type temperatures are often in Kelvin; convert so "+2C" scenarios are interpretable
for (v in intersect(c("tmax", "tmean", "dew", "utci"), names(wx))) {
  if (median(wx[[v]], na.rm = TRUE) > 100) {
    message(v, " looks like Kelvin -> converting to Celsius")
    wx[[v]] <- wx[[v]] - 273.15
  }
}
cat("Weather summary (tmax, tmean):\n"); print(summary(wx[, c("tmax", "tmean")]))

# Lags / rolling means computed on the CONTIGUOUS weather series (before joining to the
# demand table, which has missing days), so "yesterday" really means yesterday.
wx <- wx %>%
  arrange(date) %>%
  mutate(tmax_l1    = lag(tmax, 1),
         tmax_roll3 = as.numeric(stats::filter(tmax, rep(1/3, 3), sides = 1)))
saveRDS(wx, "data/processed/wx_daily.rds")   # full 1979-2024 weather for scenario resampling

# ---- merge, clean, derive ---------------------------------------------------------
tn <- demand %>%
  left_join(wx, by = "date") %>%
  arrange(date) %>%
  mutate(
    energy_MU = ifelse(energy_MU > 0, energy_MU, NA_real_),
    peak_met  = ifelse(peak_met  > 0, peak_met,  NA_real_),
    unmet_use = if (UNMET_NA_AS_ZERO) replace_na(peak_unmet, 0) else peak_unmet,
    peak_demand = peak_met + unmet_use           # underlying demand proxy
  ) %>%
  select(-unmet_use) %>%
  add_calendar()

# ---- data-quality report --------------------------------------------------------
stopifnot(!anyDuplicated(tn$date))
full_seq <- seq(min(tn$date), max(tn$date), by = "day")
cat("Date range:", format(min(tn$date)), "to", format(max(tn$date)),
    "| missing calendar days:", length(setdiff(full_seq, tn$date)), "\n")
cat("NA counts:\n"); print(colSums(is.na(tn[, c("energy_MU", "peak_met", "peak_unmet",
                                                "peak_demand", names(WX_PATTERNS))])))
cat("Share of days with unmet peak > 0:",
    round(mean(tn$peak_unmet > 0, na.rm = TRUE), 3), "\n")

saveRDS(tn, "data/processed/tn_daily.rds")
write_csv(tn, "data/processed/tn_daily.csv")
cat("Saved data/processed/tn_daily.rds (", nrow(tn), "rows )\n")