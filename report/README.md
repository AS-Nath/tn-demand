[← Back to the project README](../README.md)

# Quantifying Uncertainty in Electricity Demand: A Monte Carlo Study of Tamil Nadu’s Power System

## 1. Question

A point forecast such as “peak demand will be 18,000 MW” hides the quantity a planner actually needs: the chance that demand exceeds available capacity. Let $D_{\text{peak}}$ be the annual maximum of daily peak demand and $C$ a threshold taken from the data (a previous record, a high percentile, or the outcome that actually occurred). We estimate


$$
p(C) = P(D_{\text{peak}} > C)
$$


by simulating many plausible daily paths from a weather-conditioned regression, and we ask whether the resulting probability statements are calibrated.

Two modelling choices are treated as part of the question, not as footnotes:

1. **Supplied energy is not automatically underlying demand.** The files distinguish peak power met from peak unmet. We analyse them separately.
2. **Monte Carlo error is not the same as statistical uncertainty.** We report a binomial sampling SE, the spread of risk across parameter draws, and the spread across model specifications.

## 2. Data

Daily Tamil Nadu series, 1 January 2013 – 28 April 2024:

| Series | Units | Role |
|---|---|---|
| Corrected daily energy met | MU (GWh) | Energy supplied |
| Corrected peak met | MW | Peak power supplied |
| Peak unmet | MW | Peak deficit vs demand |
| ERA5-style weather | °C, etc. | `tmax`, `tmean`, dewpoint, solar, cloud, UTCI |

Temperatures in the state file are Kelvin and were converted to Celsius. Conflicting duplicate scrapes on a handful of dates were set to NA. About 236 calendar days are missing from the demand files (largest gap: January–March 2013).

**Underlying peak** is defined as `peak_demand = peak_met + peak_unmet` when unmet is observed. After 2015 this distinction almost vanishes:

| Year | Share of days with unmet $> 0$ | Max unmet (MW) |
|---|---|---|
| 2013 | 99.6% | 2,700 |
| 2014 | 68.1% | 2,910 |
| 2015 | 45.6% | 2,333 |
| 2016 | 0.3% | 30 |
| 2017–2024 | ≈ 0 | 0 (one 671 MW day in 2018) |

On days with unmet $> 0$, the deficit is typically about 6% of underlying peak (mean 7.1%, max 23%). A logistic model for $P(\text{unmet}>0)$ is dominated by a time trend ($z = -17.6$): the grid’s ability to meet peak, not temperature, is the main change. **Decision:** fit demand models from 2016-01-01 onward, where `peak_met` ≈ `peak_demand`. Using 2013–2015 “peak met” as if it were demand would censor the upper tail.

2024 is only a partial year (through 28 April), so the held-out evaluation year is **2023**.

## 3. Historical structure (Phase 1)

- **Level and trend.** Mean daily peak rose from about 13.4 GW (2016) to 16.1 GW (2023). Annual maxima increased by an estimated 386 MW per year (95% CI 155–617) on years with $\ge 300$ observations.
- **Season.** Kruskal–Wallis test of peak demand by month: $\chi^2_{11} = 483$, $p < 10^{-16}$. The upper tail is concentrated in February–April (April p95 ≈ 19.4 GW). November is the quietest.
- **Calendar.** ANOVA on $\log(\text{peak})$: weekday and holiday effects are both highly significant. Sundays and public holidays run well below weekdays.
- **Weather.** Spearman $\rho(\text{peak}, t_{\max}) = 0.31$ ($p < 10^{-16}$). Solar radiation is similarly associated ($\rho = 0.33$); cloud cover is negatively associated ($\rho = -0.21$).
- **Persistence.** Raw demand is strongly autocorrelated. This matters for simulating *annual maxima*: clustered high days are not 365 independent draws.

The ten highest days in the whole file are all in **April 2024** (max 20,583 MW on 26 April), already above 2023’s maximum. Any 2023-trained model that ignores growth will be too low for 2024.

## 4. Regression model (Phase 2)

All models are for $Y_t = \log(D_t)$, where $D_t$ is daily peak demand, trained on 2016-01-01 to 2021-12-31 ($n=2070$) and tested on 2022–2024 YTD ($n=841$).


$$
Y_t = \beta_0 + \beta_1 t_{\text{years}} + \text{dow}_t + \text{holiday}_t + \text{festival}_t + \text{covid}_t + \text{lockdown}_t
+ \sum_{k=1}^{3} (a_k \sin k\omega t + b_k \cos k\omega t) + f(\text{weather}_t) + \varepsilon_t.
$$


| Model | Weather | Adj. $R^2$ | Test RMSE (MW) | Test MAPE |
|---|---|---|---|---|
| M1 calendar | none | 0.586 | 1,449 | 7.16% |
| M2 linear | $t_{\max}$ | 0.692 | 1,275 | 6.36% |
| M3 spline | $ns(t_{\max},4)+$ dew | 0.700 | 1,248 | 6.24% |
| **M4 full (chosen)** | spline $t_{\max}$, 3-day roll, dew, solar, cloud | **0.717** | **1,202** | **6.01%** |
| M5 GAM | smooth $t_{\max}$, dew | 0.705 | 1,234 | 6.17% |

Newey–West HAC standard errors (14 lags) on M4:

- Coefficients are on the log scale; effects below are $100(e^{\hat\beta}-1)$ percent. Trend: +1.96% per year ($t=12.5$).
- Sunday -9.2%, fixed-date public holiday -9.5%, Diwali window -16.3%, and the hard 2020 lockdown -23.1% (lockdown indicator plus the broader COVID indicator).
- Natural spline in $t_{\max}$ and the 3-day rolling maximum are all positive and significant; dewpoint is negative (humid, cloudy days after controlling for temperature).

**Diagnostics (training residuals).** Skewness $-0.60$, kurtosis $5.75$ (heavier tails than Gaussian). AR(1) $= 0.67$. Breusch–Pagan and Durbin–Watson / Ljung–Box tests all reject iid homoscedastic errors ($p \ll 0.001$). Residual SD is largest in May (0.070 on the log scale) and smallest in February (0.033). The Monte Carlo engine therefore allows month-specific scale, AR(1) on standardised residuals, and Gaussian / $t_5$ / residual-bootstrap shocks.

The log model plus simulated $\varepsilon$ is **not** the same as $\exp(\hat\mu)$. The point forecast $\exp(\hat\mu)$ is a conditional median under symmetric errors; exceedance probabilities use the full simulated distribution.

## 5. Monte Carlo (Phase 3)

**Design.** Fit M4 on 2016–2022. For 2023, take the observed (or scenario) weather and calendar. For parameter-uncertainty runs, residual *block* bootstrap (block length 14) the training series, refit, then simulate $M$ error paths. Statistic: share of paths whose **maximum** exceeds $C$. Thresholds are data-anchored: training 99th percentile (16,668 MW), 2022 record (17,563 MW), 2023 observed max (19,387 MW).

### 5.1 Exceedance curve, 2023, observed weather

With parameter uncertainty, $N = B \times M = 50{,}000$ paths:

| $C$ (MW) | $\hat p$ | MC SE (binomial) | Parameter SD | 95% interval for $\hat p$ |
|---|---|---|---|---|
| 17,563 (2022 record) | 0.999 | 0.0002 | 0.004 | 0.998–0.999 |
| 19,000 | 0.430 | 0.0022 | 0.134 | 0.420–0.439 |
| 19,387 (what happened) | 0.232 | 0.0019 | 0.105 | 0.225–0.239 |
| 20,000 | 0.073 | 0.0012 | 0.047 | 0.069–0.076 |

Breaking the 2022 record in 2023 was, given the fitted trend and 2023’s weather, essentially certain. The realised 2023 maximum sat around the 77th percentile of the simulated annual-max distribution ($\hat p = 0.23$ of exceeding it), i.e. a high but not shocking outcome *under this model*.

**Point forecast versus distribution.** Reading the exceedance curve backwards (linear interpolation on the 250 MW grid), the simulated 2023 annual maximum has median about 18.9 GW and a 90% interval of roughly 18.0–20.2 GW. The maximum of the deterministic path $\exp(\hat\mu_t)$ is only about 17.1 GW (dashed line in Fig. 10): below even the 5th percentile of the simulated annual maximum. The maximum of a noisy path is systematically higher than the maximum of its mean path, so a single-number forecast built from $\hat\mu_t$ understates peak risk by about 1.8 GW at the median, before any growth surprise.

Parameter SD is an order of magnitude larger than binomial MC error at these thresholds: extra simulations would not change the story; extra data or a better trend model would.

### 5.2 Weather sensitivity

Additive shifts of 2023’s observed $t_{\max}$ (and lagged/rolling temperature), holding other covariates as observed:

| Scenario | $P(\max > 19{,}387)$ |
|---|---|
| Observed 2023 weather | 0.23 |
| +1 °C | 0.41 |
| +2 °C | 0.63 |
| +3 °C | 0.81 |
| Random historical year (climatology) | 0.43 |
| Climatology +2 °C | 0.81 |

A two-degree warming of the same calendar year more than doubles the estimated probability of exceeding 2023’s actual peak. Replacing 2023 weather with a random historical year *raises* the probability (0.23 → 0.43): 2023’s weather, given the model, was not the hottest plausible draw.

### 5.3 Model uncertainty

At $C = 19{,}387$ MW, plug-in $\hat p$ across specification × error law × AR(1) ranges from **0.16** (M2, bootstrap, AR) to **0.56** (M3, $t_5$, iid). Weather specification (M2 vs M4) moves the estimate only a few points. The error law and whether residuals are iid vs AR(1) move it by a factor of two to three:

- iid paths overstate annual-max risk (365 nearly independent chances at the tail);
- $t_5$ shocks fatten the tail relative to Gaussian;
- residual bootstrap *thins* the upper tail (training residuals are left-skewed).

### 5.4 Monte Carlo error check

At $C = 19{,}500$ MW ($\hat p \approx 0.20$), the SD of $\hat p$ over 30 independent reruns of $N=1000$ plug-in paths was **0.0107**, against the binomial formula $\sqrt{\hat p(1-\hat p)/N} = 0.0124$. Sampling error is behaving as theory says.

### 5.5 Seasonal risk

Let $C = 19{,}000$ MW, the level the 2023 annual maximum exceeds about half the time (plug-in curve). For a 2023-like year (trained through 2022), $P(\text{at least one day above } C)$ by month is **0.31 in May** (90% band 0.29–0.32), **0.15 in April**, 0.04 in June, 0.005 in March and July, and below 0.005 in every other month. The parameter-uncertainty spread is largest in May and April (0.12 and 0.09). High-demand risk is a **late-spring** phenomenon in this climate, concentrated in April and May. May ranks above April because of May’s higher residual variability and the 2023 weather. This replaces an earlier version that used the training 95th percentile of daily demand as the threshold, which gave probabilities near 1 in most months and said nothing about seasonality.

### 5.6 Forward look (next 365 days after 28 April 2024)

Refit on all modelling data; weather resampled from complete historical years (1979–2023). Let $C_{\text{all}} = 20{,}583$ MW, the true record (26 April 2024). With climatological weather, $P(\max > C_{\text{all}}) \approx 0.73$ (parameter-uncertainty SD 0.19); with climatology +2 °C, $\approx 0.99$. Exceeding 1.05 times the record (21,612 MW) has probability 0.22 under climatology and 0.66 with +2 °C. Under the linear trend, new records are the default, not the surprise. These probabilities depend on the trend continuing, which Section 6 suggests is likely too low in a high-growth year.

### 5.7 Comparison with a Gumbel model for annual maxima

A stationary Gumbel fit to seven annual maxima (2016–2022) gives $P(\max_{2023} > 19{,}387) \approx 0.002$. Allowing a linear location trend (slope ≈ 286 MW/year, better AIC) still gives only 0.003. Seven yearly points cannot see 2023’s weather or the daily structure that the regression uses. The Gumbel answers a different question (“how unusual is this annual max in a short list of annual maxes?”) and **understates** growth-plus-weather risk. That contrast is the reason to do the Monte Carlo.

## 6. Are the probabilities reliable? (Phase 4)

Rolling origin: for each year $Y = 2019, \ldots, 2023$, fit on data before 1 January $Y$, simulate that year’s days given **observed** weather, and check whether a nominal $1-\alpha$ central predictive interval contains the observation. Five error models; 1,754 test days.

**Pooled over all five years, a nominal 95% interval contains only about 80% of days** (0.79–0.80 across variants; 0.64–0.68 in March–May). Kupiec tests reject a 5% miss rate ($p \approx 0$), Christoffersen tests reject independence of misses, and the PIT histograms are U-shaped. Taken at face value this says the intervals are far too narrow. The pooled figure, however, mixes very different years. Coverage by year for the headline variant (AR(1), month-specific SD, parameter uncertainty):

| Year | Days | Coverage of nominal 95% | Actual below 2.5% quantile | Actual above 97.5% quantile |
|---|---|---|---|---|
| 2019 | 361 | 0.90 | 9.4% of days | 0.3% |
| 2020 | 316 | **0.51** | **48.7%** | 0.3% |
| 2021 | 352 | **0.76** | 5.1% | **18.8%** |
| 2022 | 361 | 0.92 | 1.4% | 6.1% |
| 2023 | 364 | 0.85 | 0.5% | **14.3%** |

Coverage of a nominal 95% interval for all five variants, by subset of years:

| Variant | All years | Excl. 2020–21 | 2022 | 2023 |
|---|---|---|---|---|
| iid normal, plug-in | 0.80 | 0.91 | 0.95 | 0.87 |
| AR(1) normal, plug-in | 0.80 | 0.91 | 0.95 | 0.87 |
| AR(1) + month-specific SD, plug-in | 0.79 | 0.89 | 0.91 | 0.84 |
| + parameter uncertainty | 0.80 | 0.89 | 0.92 | 0.85 |
| + residual bootstrap | 0.79 | 0.88 | 0.91 | 0.81 |

For the headline variant, March–May coverage is 0.67 pooled but 0.92 once 2020–21 are excluded (0.98 in 2022, 0.82 in 2023). Excluding 2020–21 (1,086 days), the 50/80/90/95/99% intervals cover 0.39/0.69/0.81/0.89/0.97 of days: still too narrow, by six to eleven points at the 50–95% levels (against fourteen to eighteen points when pooled).

**Where the misses come from.**

- **2020 (COVID lockdown).** Demand fell *below* the 2.5% forecast quantile on 49% of days. Models trained on 2016–2019 contain no pandemic, so no weather-and-calendar regression could anticipate it. This is an unforecastable regime shock, not a failure of the error model.
- **2021.** Misses are mostly on the high side (19% of days). This is consistent with the single COVID indicator, estimated mainly from 2020, applying too large a reduction to early 2021 as demand recovered (a hypothesis; not tested separately).
- **2022.** Close to nominal (0.92).
- **2023.** Every excess miss is on the high side: actual demand was above the 97.5% quantile on 14% of days. This is the growth acceleration that a linear log-trend fitted on earlier years cannot see, and it is the failure that matters for peak-risk statements.

Other observations: AR(1) barely moves *daily* coverage (a stationary AR has the same marginal variance as the iid model), though it matters for annual maxima (Section 5.3). The residual-bootstrap variant has the lowest 2023 coverage (0.81) because it thins the upper tail.

**Verdict.** The pooled “80% instead of 95%” should not be read as a general property of the method. In ordinary years nominal 95% intervals cover about 89–92% (modestly too narrow); in COVID-disrupted years nothing trained on earlier data is calibrated; and in the one post-COVID growth year, errors are one-sided and push upward. A planner using $p(C)$ from M4 should treat tail probabilities as likely biased *low* in a year of accelerating demand. The Monte Carlo is honest about sampling error and useful for comparing scenarios, but its absolute probabilities inherit the bias of the mean model.

## 7. Answer to the research question

1. **High-demand days are structured**, not white noise: season, weekday, holidays, temperature, and a strong trend.
2. **Unmet peak is a 2013–2015 phenomenon** in this file. From 2016, supplied peak is a reasonable demand proxy; before that it is censored.
3. **Given 2023 weather and the 2016–2022 model, exceeding the 2022 record was ~certain** ($\hat p \approx 1$); exceeding the eventual 19,387 MW outcome had probability about **0.23**. Heating the same year by 2 °C raises that to about **0.63**. A point forecast of the year’s peak (about 17.1 GW) sits below the 5th percentile of the simulated annual maximum (median 18.9 GW).
4. **Specification of the error process changes risk more than adding extra weather terms.** Ignoring autocorrelation inflates annual-max probabilities.
5. **Monte Carlo sampling error is small** next to parameter SD and model choice.
6. **Calibration is acceptable in ordinary years and fails in disrupted ones.** Pooled coverage of nominal 95% intervals is 0.80, but 0.89 excluding 2020–21 (0.92 in 2022, 0.85 in 2023). The 2020 failure is the COVID lockdown, which no earlier training window could foresee; the 2023 failure is one-sided and consistent with growth running ahead of a linear trend. Tail probabilities such as $p(C)$ are therefore likely biased low in a year of accelerating demand.

## 8. Limitations

- Linear $\log$-trend cannot capture the 2022–24 acceleration; a change-point or GLM with a more flexible trend would likely raise far-tail $p(C)$.
- COVID is represented by two indicator variables (a constant effect for March 2020–June 2021 plus an extra effect for the March–May 2020 lockdown). Rolling-origin forecasts for 2020 cannot foresee the lockdown, and 2021 forecasts inherit a possibly too-large constant adjustment.
- The highest day in the file (20,583 MW on 26 April 2024) is dropped from the modelling *data* because solar radiation is missing that day. It is used as the forward-look threshold in 5.6 but does not enter the model fit.
- Weather is exogenous. No feedback from demand to operations.
- Missing days slightly lower observed annual maxima.
- Festival dates are a ±1-day window around Diwali, plus a fixed holiday list.
- National hourly renewables (Mendeley) were not used: they are not Tamil Nadu demand.

## 9. How to reproduce

```r
source("R/run_all.R")
```

Figures `../figures/01_*.png` to `16_*.png`; tables `../outputs/*.csv`. Configuration: `R/specs.R`.

## References

Hunt, K. & Bloomfield, H. (2025). Daily electricity demand data for Indian states (2014–2024). Zenodo. https://doi.org/10.5281/zenodo.14983362

The dataset record links a companion paper, “Building and explaining data-driven energy demand models for Indian states” (see the Zenodo record for the full citation).