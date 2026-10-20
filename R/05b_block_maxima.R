# R/05b_block_maxima.R -- Gumbel / GEV comparison with the regression Monte Carlo
# Annual maxima are a natural alternative estimator of P(year-max > C).
# A linear trend in the location parameter is needed because demand is growing.
source("R/mc_functions.R")

d <- readRDS("data/processed/tn_daily.rds") %>% prep_model_data()
full_years <- d %>% count(year) %>% filter(n >= 300) %>% pull(year)
TEST_YEAR  <- max(full_years)
train  <- d %>% filter(year < TEST_YEAR)
target <- d %>% filter(year == TEST_YEAR)

ann <- train %>%
  group_by(year) %>%
  summarise(n = n(), M = max(y), .groups = "drop") %>%
  mutate(t = year - min(year))
C_rec <- max(train$y)
C_act <- max(target$y)
C_p99 <- unname(quantile(train$y, 0.99))

# Gumbel CDF: exp(-exp(-(x-mu)/sigma)); nonstationary location mu = a + b t
nll_gumbel <- function(par, M, t, trend = TRUE) {
  a <- par[1]; b <- if (trend) par[2] else 0; s <- tail(par, 1)
  if (s <= 0) return(1e12)
  mu <- a + b * t
  z  <- (M - mu) / s
  sum(z + exp(-z) + log(s))
}

fit_g <- function(trend) {
  start <- if (trend) c(a = min(ann$M), b = 200, s = sd(ann$M)) else c(a = mean(ann$M), s = sd(ann$M))
  opt <- optim(start, nll_gumbel, M = ann$M, t = ann$t, trend = trend, method = "L-BFGS-B",
               lower = if (trend) c(-Inf, -Inf, 1e-3) else c(-Inf, 1e-3),
               hessian = TRUE)
  k <- length(opt$par)
  list(par = opt$par, nll = opt$value, aic = 2 * k + 2 * opt$value,
       bic = k * log(nrow(ann)) + 2 * opt$value, trend = trend, conv = opt$convergence)
}

g_stat <- fit_g(FALSE)
g_tr   <- fit_g(TRUE)

t_star <- TEST_YEAR - min(ann$year)
mu_stat <- unname(g_stat$par[1])
s_stat  <- unname(g_stat$par[2])
mu_tr   <- unname(g_tr$par[1] + g_tr$par[2] * t_star)
s_tr    <- unname(g_tr$par[3])

pgumbel_ex <- function(C, mu, s) 1 - exp(-exp(-(C - mu) / s))

C_grid <- sort(unique(c(C_p99, C_rec, C_act, seq(16000, 22000, by = 250))))
tabG <- tibble(
  C = C_grid,
  p_gumbel_stationary   = pgumbel_ex(C_grid, mu_stat, s_stat),
  p_gumbel_trend        = pgumbel_ex(C_grid, mu_tr, s_tr)
)

# Overlay the Monte Carlo exceedance curve if available
if (file.exists("outputs/A_exceedance_curve.csv")) {
  mc <- read_csv("outputs/A_exceedance_curve.csv", show_col_types = FALSE) %>%
    filter(version == "with parameter uncertainty") %>%
    transmute(C, p_mc = p_hat)
  tabG <- tabG %>% left_join(mc, by = "C")
}

write_csv(tabG, "outputs/G_gumbel_vs_mc.csv")

cat("Annual maxima used for Gumbel fit:\n"); print(ann)
cat("\nStationary Gumbel  AIC", round(g_stat$aic, 1), "  mu", round(mu_stat), "sigma", round(s_stat), "\n")
cat("Trend Gumbel       AIC", round(g_tr$aic, 1),
    "  mu(2023)", round(mu_tr), "sigma", round(s_tr),
    "  location slope", round(g_tr$par[2]), "MW/year\n")
cat("\nP(max > training record", round(C_rec), "MW):\n")
cat("  stationary Gumbel:", round(pgumbel_ex(C_rec, mu_stat, s_stat), 3), "\n")
cat("  trend Gumbel:     ", round(pgumbel_ex(C_rec, mu_tr, s_tr), 3), "\n")
cat("P(max > actual", TEST_YEAR, "max", round(C_act), "MW):\n")
cat("  stationary Gumbel:", round(pgumbel_ex(C_act, mu_stat, s_stat), 3), "\n")
cat("  trend Gumbel:     ", round(pgumbel_ex(C_act, mu_tr, s_tr), 3), "\n")

long <- tabG %>%
  pivot_longer(-C, names_to = "method", values_to = "p") %>%
  filter(!is.na(p)) %>%
  mutate(method = recode(method,
                         p_gumbel_stationary = "Gumbel, stationary",
                         p_gumbel_trend = "Gumbel, linear location trend",
                         p_mc = "Regression Monte Carlo"))
ggsave("figures/16_gumbel_vs_mc.png",
  ggplot(long, aes(C, p, colour = method)) + geom_line() +
    geom_vline(xintercept = C_act, linetype = 3) +
    labs(title = paste("P(annual maximum > C) in", TEST_YEAR),
         subtitle = "block-maxima Gumbel vs weather-conditioned Monte Carlo",
         x = "C (MW)", y = "probability", colour = NULL),
  width = 9, height = 5, dpi = 150)
