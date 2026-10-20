# R/specs.R -- shared configuration, calendar features, model specifications
# Every script starts with source("R/specs.R") and is run from the project root.

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(splines)
})

# ---- CONFIG (edit here, nowhere else) ---------------------------------------
TARGET      <- "peak_demand"           # "peak_demand" (= met + unmet) or "energy_MU"
MODEL_START <- as.Date("2016-01-01")   # unmet-demand regime ends ~2016; model from here
TRAIN_END   <- as.Date("2021-12-31")   # train: MODEL_START..TRAIN_END ; test: afterwards
HOT_MONTHS  <- 3:5                     # Mar–May: highest monthly means and upper tail
CHOSEN      <- "M4_full"               # model carried into Monte Carlo

# Weather columns the models/simulations may use (created in 01_prepare_data.R)
WX_COLS   <- c("tmax", "tmean", "dew", "solar", "cloud", "tmax_l1", "tmax_roll3")
# Columns shifted in the "hotter than usual" scenarios
TEMP_COLS <- c("tmax", "tmean", "tmax_l1", "tmax_roll3")

# ---- Festival days -----------------------------------------------------------
# Diwali/Deepavali dates, written from memory -- CHECK against the Tamil Nadu
# government holiday list. The +/- 1 day window absorbs a one-day error.
.diwali <- as.Date(c("2016-10-30", "2017-10-18", "2018-11-06", "2019-10-27",
                     "2020-11-14", "2021-11-04", "2022-10-24", "2023-11-12",
                     "2024-10-31", "2025-10-20"))
FESTIVAL_DAYS <- sort(unique(c(.diwali - 1, .diwali, .diwali + 1)))

# ---- Calendar features -------------------------------------------------------
add_calendar <- function(df) {
  df %>%
    mutate(
      year    = year(date),
      month   = month(date),
      doy     = yday(date),
      md      = format(date, "%m-%d"),
      t_years = as.numeric(date - as.Date("2014-01-01")) / 365.25,
      dow     = factor(wday(date, week_start = 1), levels = 1:7,
                       labels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")),
      sin1 = sin(2 * pi * doy / 365.25), cos1 = cos(2 * pi * doy / 365.25),
      sin2 = sin(4 * pi * doy / 365.25), cos2 = cos(4 * pi * doy / 365.25),
      sin3 = sin(6 * pi * doy / 365.25), cos3 = cos(6 * pi * doy / 365.25),
      # Fixed-date holidays: New Year, Pongal window, Republic Day, Tamil New Year,
      # May Day, Independence Day, Gandhi Jayanti, Christmas.
      holiday = as.integer(md %in% c("01-01", "01-14", "01-15", "01-16", "01-17",
                                     "01-26", "04-14", "05-01", "08-15", "10-02",
                                     "12-25")),
      festival = as.integer(date %in% FESTIVAL_DAYS),
      # Approximate COVID-disruption window (lockdown + second wave) ...
      covid    = as.integer(date >= as.Date("2020-03-25") & date <= as.Date("2021-06-30")),
      # ... plus an extra effect for the hard national lockdown
      lockdown = as.integer(date >= as.Date("2020-03-25") & date <= as.Date("2020-05-31"))
    )
}

prep_model_data <- function(df) {
  df %>%
    mutate(y = .data[[TARGET]], log_y = log(y)) %>%
    filter(date >= MODEL_START,
           if_all(all_of(c("log_y", WX_COLS)), is.finite)) %>%
    arrange(date)
}

# ---- Candidate models (all on log scale) ------------------------------------
cal <- "t_years + dow + holiday + festival + covid + lockdown + sin1 + cos1 + sin2 + cos2 + sin3 + cos3"
mk  <- function(extra) as.formula(paste("log_y ~", cal, if (nzchar(extra)) paste("+", extra)))

SPECS <- list(
  M1_calendar = mk(""),
  M2_linear   = mk("tmax"),
  M3_spline   = mk("ns(tmax, 4) + dew"),
  M4_full     = mk("ns(tmax, 4) + tmax_roll3 + dew + solar + cloud")
)

# ---- Metrics -----------------------------------------------------------------
rmse <- function(a, p) sqrt(mean((a - p)^2))
mae  <- function(a, p) mean(abs(a - p))
mape <- function(a, p) 100 * mean(abs(a - p) / a)

dir.create("figures", showWarnings = FALSE)
dir.create("outputs", showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)