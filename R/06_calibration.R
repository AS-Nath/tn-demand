# R/06_calibration.R -- Phase 4: do nominal intervals achieve their nominal coverage?
# Rolling origin: for each year Y, fit on data before Y, simulate Y given its OBSERVED weather,
# compare simulated predictive distribution with what actually happened.
source("R/mc_functions.R")

d <- readRDS("data/processed/tn_daily.rds") %>% prep_model_data()
full_years <- d %>% count(year) %>% filter(n >= 300) %>% pull(year)
eval_years <- full_years[full_years >= 2019]
SPEC <- SPECS[[CHOSEN]]
levels_nom <- c(0.50, 0.80, 0.90, 0.95, 0.99)

variants <- tribble(
  ~name,                                    ~param_unc, ~ar,   ~hetero, ~error,
  "iid normal, plug-in",                    FALSE,      FALSE, FALSE,   "normal",
  "AR(1) normal, plug-in",                  FALSE,      TRUE,  FALSE,   "normal",
  "AR(1)+month sd, plug-in",                FALSE,      TRUE,  TRUE,    "normal",
  "AR(1)+month sd, +param uncertainty",     TRUE,       TRUE,  TRUE,    "normal",
  "AR(1)+month sd, +param unc, bootstrap",  TRUE,       TRUE,  TRUE,    "bootstrap")

one_year <- function(yr, v) {
  train  <- d %>% filter(date < as.Date(paste0(yr, "-01-01")))
  target <- d %>% filter(year == yr)
  r <- run_mc(SPEC, train, target, C = max(train$y), B = 120, M = 20,
              param_unc = v$param_unc, ar = v$ar, hetero = v$hetero, error = v$error,
              keep_paths = TRUE, seed = yr)
  P <- r$paths; obs <- target$y
  out <- target %>% transmute(date, month, obs = y, variant = v$name)
  out$pit <- colMeans(sweep(P, 2, obs, "<="))                 # P(sim <= obs)
  for (l in levels_nom) {
    lo <- apply(P, 2, quantile, probs = (1 - l) / 2); hi <- apply(P, 2, quantile, probs = 1 - (1 - l) / 2)
    out[[paste0("miss_", l)]] <- as.integer(obs < lo | obs > hi)
  }
  out
}

results <- map_dfr(seq_len(nrow(variants)), function(i) {
  message("variant: ", variants$name[i])
  map_dfr(eval_years, one_year, v = variants[i, ])
}) %>% arrange(variant, date)
saveRDS(results, "outputs/calibration_raw.rds")

# ---- coverage tables ------------------------------------------------------------------------------
cov_tab <- results %>%
  pivot_longer(starts_with("miss_"), names_to = "level", values_to = "miss") %>%
  mutate(level = as.numeric(sub("miss_", "", level)), hot = month %in% HOT_MONTHS) %>%
  group_by(variant, level) %>%
  summarise(coverage_all = 1 - mean(miss), coverage_hot = 1 - mean(miss[hot]), n = n(), .groups = "drop")
write_csv(cov_tab, "outputs/calibration_coverage.csv")
print(cov_tab %>% filter(level == 0.95), width = Inf)

ggsave("figures/14_reliability.png",
  ggplot(cov_tab, aes(level, coverage_all, colour = variant)) + geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_line() + geom_point() + coord_equal() +
    labs(title = "Reliability: nominal vs empirical coverage (rolling origin)", x = "nominal", y = "empirical"),
  width = 8, height = 6, dpi = 150)

# ---- formal tests: Kupiec (coverage) and Christoffersen (independence of misses) -------------------------------
xlogy <- function(x, y) ifelse(x == 0, 0, x * log(y))
kupiec <- function(miss, p) {
  n <- length(miss); x <- sum(miss); ph <- x / n
  ll <- function(q) xlogy(n - x, 1 - q) + xlogy(x, q)
  LR <- -2 * (ll(p) - ll(ph)); c(LR = LR, p_value = pchisq(LR, 1, lower.tail = FALSE))
}
christoffersen <- function(miss) {
  a <- miss[-length(miss)]; b <- miss[-1]
  n00 <- sum(a == 0 & b == 0); n01 <- sum(a == 0 & b == 1); n10 <- sum(a == 1 & b == 0); n11 <- sum(a == 1 & b == 1)
  p01 <- n01 / max(n00 + n01, 1); p11 <- n11 / max(n10 + n11, 1); p <- (n01 + n11) / (n00 + n01 + n10 + n11)
  L0 <- xlogy(n00 + n10, 1 - p) + xlogy(n01 + n11, p)
  L1 <- xlogy(n00, 1 - p01) + xlogy(n01, p01) + xlogy(n10, 1 - p11) + xlogy(n11, p11)
  LR <- -2 * (L0 - L1); c(LR = LR, p_value = pchisq(LR, 1, lower.tail = FALSE))
}
tests <- results %>% group_by(variant) %>%
  summarise(miss_rate_95 = mean(miss_0.95),
            kupiec_p = kupiec(miss_0.95, 0.05)[["p_value"]],
            christoffersen_p = christoffersen(miss_0.95)[["p_value"]], .groups = "drop")
write_csv(tests, "outputs/calibration_tests.csv"); print(tests, width = Inf)

# ---- PIT histograms: flat = calibrated; U-shape = intervals too narrow; hump = too wide -----------------------------------------
ggsave("figures/15_pit.png",
  ggplot(results, aes(pit)) + geom_histogram(bins = 10, boundary = 0, fill = "grey60", colour = "white") +
    facet_wrap(~ variant) + labs(title = "PIT histograms", x = "PIT", y = "days"),
  width = 10, height = 6, dpi = 150)
