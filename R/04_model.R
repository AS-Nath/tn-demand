# R/04_model.R -- Phase 2: regression models, diagnostics, out-of-sample comparison
# Needs: lmtest, sandwich, mgcv (mgcv ships with R)
source("R/specs.R"); library(lmtest); library(sandwich); library(mgcv)

d     <- readRDS("data/processed/tn_daily.rds") %>% prep_model_data()
train <- d %>% filter(date <= TRAIN_END)
test  <- d %>% filter(date >  TRAIN_END)
cat("Train:", nrow(train), "days | Test:", nrow(test), "days\n")

fits <- lapply(SPECS, function(f) lm(f, data = train))

eval_one <- function(nm, fit) {
  p_tr <- exp(fitted(fit)); p_te <- exp(predict(fit, newdata = test))   # back-transformed (median) forecasts
  tibble(model = nm, k = length(coef(fit)), AIC = AIC(fit), BIC = BIC(fit),
         adjR2 = summary(fit)$adj.r.squared,
         rmse_train = rmse(train$y, p_tr), rmse_test = rmse(test$y, p_te),
         mae_test = mae(test$y, p_te), mape_test = mape(test$y, p_te))
}
comp <- bind_rows(imap(fits, ~ eval_one(.y, .x)))

# GAM comparator (smooth weather effects); reported for comparison only
g1 <- gam(log_y ~ t_years + dow + holiday + festival + covid + lockdown +
            sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + s(tmax) + s(dew),
          data = train, method = "REML")
pg <- exp(predict(g1, newdata = test))
comp <- bind_rows(comp, tibble(model = "M5_gam", k = sum(g1$edf), AIC = AIC(g1), BIC = BIC(g1),
          adjR2 = summary(g1)$r.sq, rmse_train = rmse(train$y, exp(fitted(g1))),
          rmse_test = rmse(test$y, pg), mae_test = mae(test$y, pg), mape_test = mape(test$y, pg)))
print(comp, width = Inf); write_csv(comp, "outputs/model_comparison.csv")

# ---- Diagnostics for the chosen model ---------------------------------------------------------------
fit <- fits[[CHOSEN]]; res <- residuals(fit)
skew <- mean((res - mean(res))^3) / sd(res)^3
kurt <- mean((res - mean(res))^4) / sd(res)^4
cat("\nResidual skewness:", round(skew, 3), "| kurtosis:", round(kurt, 3), "(normal: 0 / 3)\n")
cat("AR(1) coefficient of residuals:", round(cor(res[-1], res[-length(res)]), 3), "\n")
print(bptest(fit))                                  # heteroscedasticity
print(dwtest(fit))                                  # lag-1 autocorrelation
print(Box.test(res, lag = 14, type = "Ljung-Box"))  # up to two weeks
cat("\nResidual SD by month (log scale):\n"); print(round(tapply(res, train$month, sd), 4))

png("figures/08_diagnostics.png", 1400, 1000, res = 130); par(mfrow = c(2, 2))
plot(fitted(fit), res, pch = 16, cex = .4, main = "Residuals vs fitted", xlab = "fitted (log MW)", ylab = "residual"); abline(h = 0, col = 2)
qqnorm(res, pch = 16, cex = .4); qqline(res, col = 2)
acf(res, lag.max = 60, main = "ACF of residuals")
plot(train$date, res, type = "l", main = "Residuals over time", xlab = "", ylab = ""); abline(h = 0, col = 2)
dev.off()

cat("\nCoefficients with Newey-West (HAC) standard errors:\n")
print(coeftest(fit, vcov = NeweyWest(fit, lag = 14, prewhite = FALSE)))

# Held-out forecast plot
test$pred <- exp(predict(fit, newdata = test))
ggsave("figures/09_test_forecast.png",
       ggplot(test, aes(date)) + geom_line(aes(y = y), colour = "grey40", linewidth = .3) +
         geom_line(aes(y = pred), colour = "firebrick", linewidth = .4) +
         labs(title = paste("Held-out period:", CHOSEN), y = "MW", x = NULL),
       width = 9, height = 4.5, dpi = 150)
