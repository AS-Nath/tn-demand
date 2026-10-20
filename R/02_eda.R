# R/02_eda.R -- Phase 1: understand historical demand
source("R/specs.R")
d <- readRDS("data/processed/tn_daily.rds")
save_plot <- function(p, name, w = 9, h = 5) ggsave(file.path("figures", name), p, width = w, height = h, dpi = 150)

# 1. Time series with 30-day rolling mean -----------------------------------------
d <- d %>% mutate(roll30 = as.numeric(stats::filter(peak_demand, rep(1/30, 30), sides = 2)))
save_plot(
  ggplot(d, aes(date)) + geom_line(aes(y = peak_demand), colour = "grey70", linewidth = .3) +
    geom_line(aes(y = roll30), colour = "firebrick") +
    labs(title = "Tamil Nadu daily peak demand (met + unmet)", y = "MW", x = NULL),
  "01_timeseries.png")

# 2. Seasonality ----------------------------------------------------------------------
save_plot(
  ggplot(d, aes(factor(month), peak_demand)) + geom_boxplot(outlier.size = .6) +
    labs(title = "Peak demand by month", x = "Month", y = "MW"), "02_monthly_box.png")
cat("\nMonthly mean / 95th percentile of peak demand:\n")
print(d %>% group_by(month) %>%
        summarise(mean = mean(peak_demand, na.rm = TRUE),
                  p95 = quantile(peak_demand, .95, na.rm = TRUE),
                  max = max(peak_demand, na.rm = TRUE)))
print(kruskal.test(peak_demand ~ factor(month), data = d))   # month effect

# 3. Long-term change ---------------------------------------------------------------------
ann <- d %>% group_by(year) %>%
  summarise(n = sum(!is.na(peak_demand)), mean = mean(peak_demand, na.rm = TRUE),
            annual_max = max(peak_demand, na.rm = TRUE))
print(ann)
fit_trend <- lm(annual_max ~ year, data = filter(ann, n >= 300))
cat("\nTrend in annual maximum (MW/year):\n"); print(summary(fit_trend)$coefficients); print(confint(fit_trend))
save_plot(ggplot(ann, aes(year, annual_max)) + geom_point() + geom_line() +
            labs(title = "Annual maximum of daily peak demand", y = "MW"), "03_annual_max.png")

# 4. Weekday / holiday effects ---------------------------------------------------------------
print(summary(aov(log(peak_demand) ~ dow + holiday, data = filter(d, is.finite(peak_demand)))))

# 5. Weather association ---------------------------------------------------------------------------
save_plot(
  ggplot(d, aes(tmax, peak_demand, colour = factor(year))) + geom_point(alpha = .35, size = .8) +
    geom_smooth(aes(group = 1), method = "loess", colour = "black", se = TRUE) +
    labs(title = "Peak demand vs daily max temperature", x = "tmax", y = "MW", colour = NULL),
  "04_demand_vs_tmax.png")
vars <- intersect(c("peak_demand", "energy_MU", "tmax", "tmean", "dew", "solar", "cloud", "utci"), names(d))
cat("\nSpearman correlations:\n"); print(round(cor(d[vars], use = "pairwise.complete.obs", method = "spearman"), 2))
print(cor.test(d$peak_demand, d$tmax, method = "spearman", exact = FALSE))

# 6. Distribution & extremes --------------------------------------------------------------------------
png("figures/05_distribution.png", 1400, 600, res = 130); par(mfrow = c(1, 2))
hist(log(d$peak_demand), breaks = 50, main = "log peak demand", xlab = "")
qqnorm(log(d$peak_demand)); qqline(log(d$peak_demand)); dev.off()
cat("\nTop 10 days by peak demand:\n")
print(d %>% arrange(desc(peak_demand)) %>% select(date, peak_met, peak_unmet, peak_demand, tmax) %>% head(10))

# 7. Autocorrelation (matters for how we simulate!) -------------------------------------------------------
png("figures/06_acf.png", 1400, 600, res = 130); par(mfrow = c(1, 2))
acf(d$peak_demand, na.action = na.pass, lag.max = 60, main = "ACF: raw demand")
acf(diff(log(d$peak_demand)), na.action = na.pass, lag.max = 60, main = "ACF: log-differences"); dev.off()

# 8. Unmet peak overview (full treatment in 03_unmet_demand.R) --------------------------------------------------
cat("\nUnmet peak: share of days > 0 by year\n")
print(d %>% group_by(year) %>% summarise(share_unmet = mean(peak_unmet > 0, na.rm = TRUE),
                                         max_unmet = max(peak_unmet, na.rm = TRUE)))
