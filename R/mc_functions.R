# R/mc_functions.R -- Monte Carlo engine
# Each simulated "path" = one plausible sequence of daily demand over the target period.
# Per parameter draw b = 1..B:
#   1. (if param_unc) residual BLOCK-bootstrap the training data and refit  -> beta*, phi*, sd*
#   2. choose weather (observed / shifted / resampled historical year)
#   3. simulate M error paths (AR(1), optionally month-specific sd, normal / t5 / bootstrap shocks)
#   4. record, for each threshold C, the share of the M paths whose maximum exceeds C  -> p_b[b, ]
# Splitting B x M lets us separate parameter uncertainty from Monte Carlo noise.
source("R/specs.R")

# ---- design matrix for new data (uses the TRAINING knots for ns()) ---------------------------
design_matrix <- function(fit, newdata) {
  tt <- delete.response(terms(fit))
  mf <- model.frame(tt, newdata, xlev = fit$xlevels, na.action = na.fail)
  model.matrix(tt, mf, contrasts.arg = fit$contrasts)
}

# ---- residual model: month-specific sd + AR(1) on standardised residuals ----------------------
resid_params <- function(res, month, hetero = TRUE) {
  months <- as.character(1:12)
  sd_m <- if (hetero) tapply(res, month, sd) else setNames(rep(sd(res), 12), months)
  sd_m <- sd_m[months]
  if (anyNA(sd_m)) stop("Some calendar months are missing from the training data.")
  sdv <- unname(sd_m[as.character(month)])
  z   <- res / sdv
  n   <- length(z)
  phi <- cor(z[-1], z[-n])
  list(sd_m = sd_m, phi = phi, z = z, sdv = sdv, u = z[-1] - phi * z[-n])
}

make_ingredients <- function(fit, train, hetero = TRUE) {
  stopifnot(nrow(train) == length(residuals(fit)))
  X <- model.matrix(fit)
  list(fit = fit, X = X, qr = qr(X), fitted = as.numeric(fitted(fit)),
       month = train$month, hetero = hetero,
       z0 = resid_params(residuals(fit), train$month, hetero))
}

block_resample <- function(x, L) {
  n <- length(x); nb <- ceiling(n / L)
  starts <- sample.int(n, nb, replace = TRUE)
  idx <- unlist(lapply(starts, function(s) ((s - 1 + 0:(L - 1)) %% n) + 1))[1:n]
  x[idx]
}

# Residual block bootstrap: keeps autocorrelation and the seasonal sd pattern
refit_boot <- function(ing, L = 14) {
  zs     <- block_resample(ing$z0$z, L)
  y_star <- ing$fitted + ing$z0$sdv * zs
  beta   <- qr.coef(ing$qr, y_star)
  beta[is.na(beta)] <- 0
  res    <- y_star - drop(ing$X %*% beta)
  c(list(beta = beta), resid_params(res, ing$month, ing$hetero))
}

# ---- standardised error paths (M paths x T days) -------------------------------------------------
sim_z <- function(par, T, M, ar = TRUE, error = "normal") {
  phi  <- if (ar) par$phi else 0
  s    <- sqrt(1 - phi^2)
  pool <- if (ar) par$u else par$z
  draw  <- function(n) switch(error,
    normal    = rnorm(n, 0, s),
    t5        = rt(n, 5) * s * sqrt(3 / 5),            # t5 rescaled to the same variance
    bootstrap = sample(pool, n, replace = TRUE))
  draw0 <- function(n) switch(error,
    normal    = rnorm(n),
    t5        = rt(n, 5) * sqrt(3 / 5),
    bootstrap = sample(par$z, n, replace = TRUE))
  Z <- matrix(0, M, T)
  Z[, 1] <- draw0(M)
  if (T > 1) for (t in 2:T) Z[, t] <- phi * Z[, t - 1] + draw(M)
  Z
}

# ---- weather scenarios ----------------------------------------------------------------------------------
build_pool <- function(train = NULL) {    # historical years with a complete weather record
  md365 <- format(seq(as.Date("2021-01-01"), as.Date("2021-12-31"), by = "day"), "%m-%d")
  readRDS("data/processed/wx_daily.rds") %>%
    add_calendar() %>%
    group_by(year) %>%
    filter(all(md365 %in% md), all(complete.cases(across(all_of(WX_COLS))))) %>%
    ungroup()
}

scenario_newdata <- function(target, scenario, pool) {
  nd <- target
  if (identical(scenario$type, "clim")) {       # one random historical year's weather per draw
    y  <- sample(unique(pool$year), 1)
    lk <- pool %>% filter(year == y)
    idx <- match(nd$md, lk$md)
    idx[is.na(idx)] <- match("02-28", lk$md)    # Feb 29 fallback
    for (v in WX_COLS) nd[[v]] <- lk[[v]][idx]
  }
  if (!is.null(scenario$delta) && scenario$delta != 0)
    for (v in TEMP_COLS) nd[[v]] <- nd[[v]] + scenario$delta
  nd
}

# ---- main simulator ----------------------------------------------------------------------------------------
# C: vector of thresholds (same units as TARGET). Statistic: P(max over target period > C).
run_mc <- function(spec, train, target, C, B = 1000, M = 50, param_unc = TRUE,
                   error = "normal", ar = TRUE, hetero = TRUE,
                   scenario = list(type = "observed", delta = 0),
                   L = 14, keep_paths = FALSE, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  fit  <- lm(spec, data = train)
  beta0 <- coef(fit); beta0[is.na(beta0)] <- 0
  ing  <- make_ingredients(fit, train, hetero)
  par0 <- c(list(beta = beta0), ing$z0)
  pool <- build_pool(train)
  T    <- nrow(target); mon <- as.character(target$month)
  p_b  <- matrix(NA_real_, B, length(C))
  paths <- if (keep_paths) matrix(NA_real_, B * M, T) else NULL

  for (b in seq_len(B)) {
    if (b == 1 || b %% 200 == 0) message("  MC draw ", b, "/", B)
    par <- if (param_unc) refit_boot(ing, L) else par0
    nd  <- scenario_newdata(target, scenario, pool)
    mu  <- drop(design_matrix(fit, nd) %*% par$beta)
    Z   <- sim_z(par, T, M, ar, error)
    sdn <- unname(par$sd_m[mon])
    D   <- exp(sweep(Z, 2, sdn, "*") + rep(mu, each = M))     # M x T, original units
    mx  <- D[, 1]
    if (T > 1) for (j in 2:T) mx <- pmax(mx, D[, j])
    p_b[b, ] <- vapply(C, function(cc) mean(mx > cc), numeric(1))
    if (keep_paths) paths[((b - 1) * M + 1):(b * M), ] <- D
  }
  list(p_b = p_b, C = C, B = B, M = M, paths = paths, fit = fit)
}

# ---- uncertainty decomposition ---------------------------------------------------------------------------------
# p_hat      : estimate of P(max > C), integrating over parameter uncertainty
# se_total   : standard error of p_hat (Monte Carlo noise + draw-to-draw variation, since draws are iid)
# mc_se_binom: pure binomial MC error sqrt(p(1-p)/N) -- the right SE if parameters are fixed
# param_sd   : spread of the *conditional* risk across plausible parameter values,
#              after removing the within-draw binomial noise  (the "model-estimation uncertainty")
summarise_p <- function(res) {
  p_b <- res$p_b; B <- nrow(p_b); M <- res$M
  p_hat <- colMeans(p_b)
  v     <- apply(p_b, 2, var)
  se    <- sqrt(v / B)
  tibble(C = res$C, p_hat = p_hat, se_total = se,
         mc_se_binom = sqrt(p_hat * (1 - p_hat) / (B * M)),
         param_sd = sqrt(pmax(0, v - colMeans(p_b * (1 - p_b)) / M)),
         lo = pmax(0, p_hat - 1.96 * se), hi = pmin(1, p_hat + 1.96 * se))
}
