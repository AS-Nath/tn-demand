# R/05_monte_carlo.R -- Phase 3: exceedance probabilities, scenarios, model uncertainty
# Run time: several minutes. Section A is cached in outputs/A_exceedance_curve.csv and
# outputs/A_mc_runs.rds; delete those two files to force a fresh run.
source("R/mc_functions.R")

d <- readRDS("data/processed/tn_daily.rds") %>% prep_model_data()
full_years <- d %>% count(year) %>% filter(n >= 300) %>% pull(year)
TEST_YEAR  <- max(full_years)                         # last complete year (2023; 2024 is partial)
train  <- d %>% filter(year < TEST_YEAR)              # everything before the evaluation year
target <- d %>% filter(year == TEST_YEAR)             # held-out year, observed weather
SPEC   <- SPECS[[CHOSEN]]
cat("MC train:", format(min(train$date)), "to", format(max(train$date)),
    "| target year:", TEST_YEAR, "(", nrow(target), "days )\n")

# Data-anchored thresholds (no arbitrary numbers). Substitute an official capacity figure if you have one.
C_rec <- max(train$y)                       # previous record (end of 2022)
C_p99 <- unname(quantile(train$y, 0.99))
C_act <- max(target$y)                      # what actually happened (retrospective check)
C_grid <- sort(unique(c(C_p99, C_rec, C_act, seq(16000, 22000, by = 250))))

# ---- A. Exceedance curve: plug-in vs with parameter uncertainty ----------------------------------------------
a_csv <- "outputs/A_exceedance_curve.csv"
a_rds <- "outputs/A_mc_runs.rds"
if (file.exists(a_csv)) {
  message("Loading cached exceedance curve from ", a_csv)
  tabA <- read_csv(a_csv, show_col_types = FALSE)
  A_plug <- if (file.exists(a_rds)) readRDS(a_rds)$A_plug else NULL
} else {
  A_full <- run_mc(SPEC, train, target, C_grid, B = 1000, M = 50, param_unc = TRUE,  seed = 1)
  A_plug <- run_mc(SPEC, train, target, C_grid, B = 1000, M = 50, param_unc = FALSE, seed = 1)
  saveRDS(list(A_full = A_full, A_plug = A_plug), a_rds)
  tabA <- bind_rows(mutate(summarise_p(A_full), version = "with parameter uncertainty"),
                    mutate(summarise_p(A_plug), version = "plug-in (parameters fixed)"))
  write_csv(tabA, a_csv)
}
print(tabA %>% filter(C %in% c(C_p99, C_rec, C_act)), width = Inf)

pl    <- tabA %>% filter(version == "plug-in (parameters fixed)")
C_med <- pl$C[which.min(abs(pl$p_hat - 0.5))]   # level the 2023 maximum exceeds about half the time
C_key <- c(C_med, C_act)

# Point forecast vs the simulated distribution of the annual maximum
mu_det <- exp(predict(lm(SPEC, data = train), newdata = target))     # naive point forecast
tabA_wp <- tabA %>% filter(version == "with parameter uncertainty", p_hat > 0.001, p_hat < 0.999)
q_annual <- approx(tabA_wp$p_hat, tabA_wp$C, xout = c(0.95, 0.5, 0.05))$y
cat("\nMax of point forecast:", round(max(mu_det)), "MW | actual", TEST_YEAR, "max:", round(C_act), "MW\n")
cat("Simulated annual max (with parameter uncertainty), 5th / 50th / 95th percentile:",
    paste(round(q_annual), collapse = " / "), "MW\n")
cat("C_med (threshold used in sections B, C, E):", C_med, "MW\n\n")

ggsave("figures/10_exceedance_curve.png",
       ggplot(tabA, aes(C, p_hat, colour = version, fill = version)) +
         geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .2, colour = NA) + geom_line() +
         geom_vline(xintercept = max(mu_det), linetype = 2) +
         geom_vline(xintercept = max(target$y), linetype = 3) +
         labs(title = paste("P(max daily demand in", TEST_YEAR, "> C)"),
              subtitle = "dashed: max of point forecast | dotted: actual max", x = "C (MW)", y = "probability"),
       width = 9, height = 5, dpi = 150)

# ---- B. Weather scenarios -------------------------------------------------------------------------------------
scen <- list(
  observed    = list(type = "observed", delta = 0),
  `+1C`       = list(type = "observed", delta = 1),
  `+2C`       = list(type = "observed", delta = 2),
  `+3C`       = list(type = "observed", delta = 3),
  climatology = list(type = "clim", delta = 0),
  `clim +2C`  = list(type = "clim", delta = 2))
tabB <- imap_dfr(scen, function(s, nm) {
  r <- run_mc(SPEC, train, target, C_key, B = 500, M = 40, scenario = s, seed = 2)
  mutate(summarise_p(r), scenario = nm)
})
write_csv(tabB, "outputs/B_weather_scenarios.csv"); print(tabB, width = Inf)

# ---- C. Model uncertainty: specification x error law x autocorrelation ------------------------------------------
grid <- expand_grid(spec = names(SPECS)[2:4], error = c("normal", "t5", "bootstrap"), ar = c(TRUE, FALSE))
tabC <- pmap_dfr(grid, function(spec, error, ar) {
  message("model grid: ", spec, " / ", error, " / ar=", ar)
  r <- run_mc(SPECS[[spec]], train, target, c(19000, C_act), B = 150, M = 25,
              error = error, ar = ar, param_unc = FALSE, seed = 3)
  mutate(summarise_p(r), spec = spec, error = error, ar = ar)
})
write_csv(tabC, "outputs/C_model_uncertainty.csv")
ggsave("figures/11_model_uncertainty.png",
       ggplot(filter(tabC, abs(C - C_act) < 1), aes(p_hat, interaction(spec, error), colour = ar)) +
         geom_pointrange(aes(xmin = lo, xmax = hi)) +
         labs(title = paste0("P(max ", TEST_YEAR, " demand > ", round(C_act), " MW)"),
              subtitle = "spread across specification, error law, and AR(1) residuals",
              x = "probability", y = NULL),
       width = 9, height = 6, dpi = 150)

# ---- D. Monte Carlo error: theory vs experiment -----------------------------------------------------------------------
C_mc <- pl$C[which.min(abs(pl$p_hat - 0.2))]          # threshold with p ~ 0.2 for a visible demo
N_per <- 1000                                          # 100 draws x 10 paths, parameters fixed => iid paths
reps <- map_dbl(1:30, ~ summarise_p(run_mc(SPEC, train, target, C_mc, B = 100, M = 10,
                                           param_unc = FALSE, seed = 100 + .x))$p_hat)
cat(sprintf("\nMC error check at C=%.0f: sd of p_hat over 30 reruns = %.4f | binomial formula = %.4f\n",
            C_mc, sd(reps), sqrt(mean(reps) * (1 - mean(reps)) / N_per)))
A_conv <- run_mc(SPEC, train, target, C_mc, B = 800, M = 20, param_unc = FALSE, seed = 6)
run_est <- tibble(N = (1:A_conv$B) * A_conv$M, p = cumsum(A_conv$p_b[, 1]) / (1:A_conv$B)) %>%
  mutate(se = sqrt(p * (1 - p) / N))
ggsave("figures/12_mc_convergence.png",
       ggplot(run_est, aes(N, p)) + geom_ribbon(aes(ymin = p - 1.96 * se, ymax = p + 1.96 * se), alpha = .2) +
         geom_line() + labs(title = "Running Monte Carlo estimate with 95% binomial CI", y = "p_hat"),
       width = 8, height = 4.5, dpi = 150)

# ---- E. Seasonal risk: P(any day in month m exceeds C_med) --------------------------------------------------------------
# C_med is a level the annual maximum exceeds about half the time. A low daily percentile
# (e.g. the training 95th) gives p ~ 1 in nearly every month and says nothing about seasonality.
C_q <- C_med
tabE <- map_dfr(1:12, function(m) {
  r <- run_mc(SPEC, train, filter(target, month == m), C_q, B = 300, M = 20, seed = 4)
  mutate(summarise_p(r), month = m, C_threshold = C_q)
})
write_csv(tabE, "outputs/E_seasonal_risk.csv"); print(select(tabE, month, p_hat, lo, hi), n = 12)
ggsave("figures/13_seasonal_risk.png",
       ggplot(tabE, aes(factor(month), p_hat)) + geom_col() +
         geom_errorbar(aes(ymin = lo, ymax = hi), width = .3) +
         labs(title = paste0("P(at least one day in month > ", round(C_q), " MW) in ", TEST_YEAR),
              subtitle = "threshold = level the annual maximum exceeds about half the time (plug-in curve)",
              x = "month", y = "probability"),
       width = 9, height = 4.5, dpi = 150)

# ---- F. Forward-looking: refit on ALL modelling data, next 365 days, climatological weather ------------------------------------
future <- add_calendar(tibble(date = seq(max(d$date) + 1, by = "day", length.out = 365)))
# True record, including 2024-04-26 (20,583 MW), which prep_model_data() drops (solar missing that day)
C_all  <- max(readRDS("data/processed/tn_daily.rds")$peak_demand, na.rm = TRUE)
tabF <- imap_dfr(scen[c("climatology", "clim +2C")], function(s, nm) {
  r <- run_mc(SPEC, d, future, c(0.95 * C_all, C_all, 1.05 * C_all), B = 500, M = 40, scenario = s, seed = 5)
  mutate(summarise_p(r), scenario = nm)
})
write_csv(tabF, "outputs/F_forward_risk.csv"); print(tabF, width = Inf)