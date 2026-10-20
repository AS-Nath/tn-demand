# R/03_unmet_demand.R -- supplied vs unmet peak, analysed separately
# Question: is peak_met a good proxy for underlying demand, or is censoring material?
source("R/specs.R"); library(splines)
d <- readRDS("data/processed/tn_daily.rds")
u <- d %>% filter(!is.na(peak_unmet), !is.na(peak_met)) %>%
  mutate(has_unmet = as.integer(peak_unmet > 0), unmet_ratio = peak_unmet / (peak_met + peak_unmet))

cat("Days with data:", nrow(u), "| days with unmet > 0:", sum(u$has_unmet),
    sprintf("(%.1f%%)\n", 100 * mean(u$has_unmet)))
cat("\nBy year:\n")
print(u %>% group_by(year) %>%
        summarise(share = mean(has_unmet),
                  mean_when_unmet = mean(peak_unmet[has_unmet == 1]),
                  max_unmet = max(peak_unmet),
                  max_ratio = max(unmet_ratio)))
cat("\nUnmet as share of underlying peak (summary over days with unmet > 0):\n")
print(summary(u$unmet_ratio[u$has_unmet == 1]))

# Part 1: WHEN does unmet demand occur?  (logistic regression)
uu <- u %>% filter(is.finite(tmax))
if (sum(uu$has_unmet) >= 30) {
  g <- glm(has_unmet ~ t_years + dow + covid + ns(tmax, 3) + sin1 + cos1 + sin2 + cos2,
           data = uu, family = binomial)
  print(summary(g)); print(drop1(g, test = "Chisq"))
} else message("Fewer than 30 days with unmet demand -- logistic model skipped.")

# Part 2: HOW MUCH, when it occurs?  (compare positive-part distributions)
pos <- u$peak_unmet[u$peak_unmet > 0]
if (length(pos) >= 30) {
  fl <- MASS::fitdistr(pos, "lognormal")
  fg <- tryCatch(MASS::fitdistr(pos, "gamma"), error = function(e) NULL)
  cat("\nLognormal AIC:", 2 * 2 - 2 * fl$loglik, "\n")
  if (!is.null(fg)) cat("Gamma AIC:    ", 2 * 2 - 2 * fg$loglik, "\n")
}

p <- u %>% group_by(year, month) %>% summarise(share = mean(has_unmet), .groups = "drop") %>%
  ggplot(aes(month, factor(year), fill = share)) + geom_tile() +
  scale_x_continuous(breaks = 1:12) + labs(title = "Share of days with unmet peak > 0", y = NULL)
ggsave("figures/07_unmet_heatmap.png", p, width = 8, height = 4.5, dpi = 150)

cat("\nDECISION GUIDE: if unmet days are rare and the ratio is tiny, peak_met ~ demand.\n",
    "Otherwise keep TARGET = 'peak_demand' (met + unmet) and report results for 'peak_met' as a sensitivity check.\n")
