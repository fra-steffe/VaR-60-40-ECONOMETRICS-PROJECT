#  Step 2 - Volatility models and one-day VaR (paper, Section 4).

#  Estimates AR(1)-GARCH(1,1) marginals for SPY and TLT, with Gaussian and
#  Student-t innovations, and a DCC(1,1) correlation on top. A rolling
#  1,000-day window, re-estimated every 21 days, gives one-step-ahead
#  covariance forecasts; the CCC model uses the same marginals with the
#  correlation held constant. Portfolio sigma is turned into 5% and 1% VaR
#  with Gaussian, Student-t and semiparametric (empirical) quantiles, plus a
#  historical-simulation benchmark. Also produces Table 6 and Figure 3.
#
#  Reproducibility: in a few windows of the Student-t branch (2013, 2015,
#  2024; none in 2022) the default optimizer fails and a fallback with random
#  restarts is used, so a re-run can differ slightly there. The CSVs in
#  results/ are the ones behind the paper.

library(rugarch)   # univariate GARCH models (the marginals)
library(rmgarch)   # DCC on top of the rugarch marginals
library(xts)       # date-indexed series, required by rmgarch
library(fGarch)    # standardized Student-t: stdFit (estimation), qstd (quantile)


# Parameters
INPUT_FILE <- "results/returns_spy_tlt_6040.csv"
RESULTS_DIR <- "results"
FIGURES_DIR <- "figures"

WEIGHTS <- c(0.60, 0.40)   # 60/40 portfolio: SPY, TLT
WINDOW_SIZE <- 1000        # rolling estimation window (~4 years)
REFIT_EVERY <- 21          # re-estimate every ~month, fixed parameters in between
TAIL_WINDOW <- 500         # past standardized returns used for the t and semiparametric tails
HS_WINDOW <- 250           # past raw returns used by historical simulation


#1. Data and model specifications 

# Load the Step 1 returns. rmgarch needs the two assets as a date-indexed
# (xts) series; the portfolio return is kept for the VaR step.
load_returns <- function(path) {
  returns <- read.csv(path)
  returns$Date <- as.Date(returns$Date)
  asset_returns <- xts(returns[, c("r_SPY", "r_TLT")], order.by = returns$Date)
  list(returns = returns, asset_returns = asset_returns)
}


# Marginal model of each asset: AR(1) mean, because Ljung-Box rejects on
# return levels (Table 3), and GARCH(1,1) variance, for the volatility
# clustering in squared returns. The innovation distribution is the first
# axis of the 2x2 design: "norm" (Gaussian) or "std" (Student-t).
marginal_spec <- function(distribution) {
  ugarchspec(
    mean.model = list(armaOrder = c(1, 0), include.mean = TRUE),
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    distribution.model = distribution
  )
}


# DCC(1,1) on top of the same marginal model for both assets.
dcc_spec <- function(marginal) {
  both_assets <- multispec(replicate(2, marginal))
  dccspec(uspec = both_assets, dccOrder = c(1, 1), distribution = "mvnorm")
}


# 2. Full-sample estimation: Table 6 and residual diagnostics
# Not used for the VaR: one estimation on
# the whole 2004-2024 sample, to report typical parameter values and check
# the AR(1)-GARCH(1,1) specification on its residuals.

fit_marginals_full_sample <- function(asset_returns, spec_normal, spec_t) {
  list(
    spy_normal = ugarchfit(spec_normal, asset_returns[, "r_SPY"]),
    spy_t = ugarchfit(spec_t, asset_returns[, "r_SPY"]),
    tlt_normal = ugarchfit(spec_normal, asset_returns[, "r_TLT"]),
    tlt_t = ugarchfit(spec_t, asset_returns[, "r_TLT"])
  )
}


# show() prints each fit with its weighted Ljung-Box tests. These are the
# diagnostics discussed in Section 4 of the paper.
print_marginal_diagnostics <- function(marginal_fits) {
  for (name in names(marginal_fits)) {
    cat("\n\n========== Marginal:", name, "==========\n")
    show(marginal_fits[[name]])
  }
}


# Table 6, panel A. "shape" is the Student-t degrees of freedom: it exists
# only for the t fits, so it is NA for the Gaussian ones.
marginal_parameter_table <- function(marginal_fits) {
  parameter_names <- c("mu", "ar1", "omega", "alpha1", "beta1")
  parameters <- data.frame(
    SPY_norm = coef(marginal_fits$spy_normal)[parameter_names],
    SPY_t = coef(marginal_fits$spy_t)[parameter_names],
    TLT_norm = coef(marginal_fits$tlt_normal)[parameter_names],
    TLT_t = coef(marginal_fits$tlt_t)[parameter_names]
  )
  parameters["shape", "SPY_t"] <- coef(marginal_fits$spy_t)["shape"]
  parameters["shape", "TLT_t"] <- coef(marginal_fits$tlt_t)["shape"]
  round(parameters, 6)
}


# Table 6, panel B: the DCC parameters a (reaction to joint shocks) and
# b (persistence), estimated on the whole sample for each branch.
dcc_parameter_table <- function(spec_normal, spec_t, asset_returns) {
  fit_normal <- dccfit(spec_normal, data = asset_returns)
  fit_t <- dccfit(spec_t, data = asset_returns)
  parameter_names <- c("[Joint]dcca1", "[Joint]dccb1")
  parameters <- data.frame(
    DCC_norm = coef(fit_normal)[parameter_names],
    DCC_t = coef(fit_t)[parameter_names]
  )
  rownames(parameters) <- c("dcc_a", "dcc_b")
  round(parameters, 6)
}


#3. DCC estimation on one block

# Estimate the DCC on one block and forecast the covariance matrix one day
# ahead for each of the block's last n_forecast days.
# dccfit() runs the two-stage estimation of Engle (2002): first the two GARCH
# marginals by maximum likelihood, then the correlation parameters a and b
# given the standardized residuals. out.sample keeps the last n_forecast days
# out of the estimation; dccforecast() forecasts each of them one step ahead
# with the parameters fixed, using only data up to the previous day.
fit_and_forecast <- function(spec, block_data, n_forecast, solver) {
  fit <- dccfit(spec, data = block_data, out.sample = n_forecast, solver = solver)
  forecast <- dccforecast(fit, n.ahead = 1, n.roll = n_forecast - 1)
  list(fit = fit, forecast = forecast)
}


# Some windows make the default optimizer fail, typically the Student-t
# marginal of TLT. Solvers are tried in order: "solnp", then "gosolnp"
# then "nlminb". If all three fail, the block is recorded as
# non-converged and its days are left missing.
estimate_dcc_block <- function(spec, block_data, n_forecast) {
  for (solver in c("solnp", "gosolnp", "nlminb")) {
    result <- tryCatch(fit_and_forecast(spec, block_data, n_forecast, solver),
                       error = function(e) NULL)
    if (!is.null(result)) {
      return(result)
    }
  }
  NULL   # all three solvers failed
}


#  4. From the covariance forecast to portfolio sigma (CCC and DCC)

# Constant correlation of the CCC model. Setting a = b = 0 in the DCC recursion leaves
# exactly this constant matrix, so the CCC needs no separate estimation and
# shares the DCC's marginals: comparing the two isolates the correlation.
constant_correlation <- function(fit) {
  standardized_residuals <- residuals(fit) / sigma(fit)
  correlation_matrix <- cor(standardized_residuals)
  correlation_matrix[1, 2]
}


# Portfolio volatility, sigma_p = sqrt(w' Sigma w).
portfolio_sigma <- function(covariance, weights) {
  sqrt(as.numeric(t(weights) %*% covariance %*% weights))
}


# For each forecast day of a block: the two marginal volatilities, and the
# portfolio sigma under the DCC (forecast covariance as it is) and under the
# CCC (same volatilities, correlation replaced by the constant one).
block_sigmas <- function(fit, forecast, weights) {
  covariances <- rcov(forecast)   # one forecast matrix per day, named by its date
  rho_constant <- constant_correlation(fit)
  n_days <- length(covariances)

  sigma_spy <- numeric(n_days)
  sigma_tlt <- numeric(n_days)
  sigma_p_dcc <- numeric(n_days)
  sigma_p_ccc <- numeric(n_days)

  for (day in 1:n_days) {
    covariance_dcc <- covariances[[day]][, , 1]
    sigma_spy[day] <- sqrt(covariance_dcc[1, 1])
    sigma_tlt[day] <- sqrt(covariance_dcc[2, 2])

    covariance_spy_tlt <- rho_constant * sigma_spy[day] * sigma_tlt[day]
    covariance_ccc <- matrix(c(sigma_spy[day]^2, covariance_spy_tlt,
                               covariance_spy_tlt, sigma_tlt[day]^2),
                             nrow = 2)

    sigma_p_dcc[day] <- portfolio_sigma(covariance_dcc, weights)
    sigma_p_ccc[day] <- portfolio_sigma(covariance_ccc, weights)
  }

  data.frame(date = as.Date(names(covariances)),
             sigma_spy = sigma_spy,
             sigma_tlt = sigma_tlt,
             sigma_p_dcc = sigma_p_dcc,
             sigma_p_ccc = sigma_p_ccc)
}


# 5. Rolling out-of-sample engine

# A 1,000-day estimation window moves forward 21 days at a time. At each step
# both branches (Gaussian and Student-t) are re-estimated on the same window
# and forecast the next 21 days with fixed parameters. Joined together, the
# forecasts of all blocks form a continuous out-of-sample series that never
# uses future data.
rolling_engine <- function(spec_normal, spec_t, asset_returns, weights,
                           window_size, refit_every) {
  n_total <- nrow(asset_returns)
  block_starts <- seq(1, n_total - window_size, by = refit_every)
  n_blocks <- length(block_starts)

  blocks_normal <- vector("list", n_blocks)
  blocks_t <- vector("list", n_blocks)
  failed_normal <- integer(0)
  failed_t <- integer(0)

  for (block in 1:n_blocks) {
    window_start <- block_starts[block]
    block_end <- min(window_start + window_size + refit_every - 1, n_total)
    block_data <- asset_returns[window_start:block_end, ]
    n_forecast <- nrow(block_data) - window_size

    estimate_normal <- estimate_dcc_block(spec_normal, block_data, n_forecast)
    if (is.null(estimate_normal)) {
      failed_normal <- c(failed_normal, block)
    } else {
      blocks_normal[[block]] <- block_sigmas(estimate_normal$fit,
                                             estimate_normal$forecast, weights)
    }

    estimate_t <- estimate_dcc_block(spec_t, block_data, n_forecast)
    if (is.null(estimate_t)) {
      failed_t <- c(failed_t, block)
    } else {
      blocks_t[[block]] <- block_sigmas(estimate_t$fit, estimate_t$forecast, weights)
    }

    cat("block", block, "of", n_blocks, "\r")
  }
  cat("\n")

  # Failed blocks are NULL and rbind() skips them: their days stay missing.
  list(normal = do.call(rbind, blocks_normal),
       t = do.call(rbind, blocks_t),
       failed_normal = failed_normal,
       failed_t = failed_t,
       n_blocks = n_blocks)
}


# Which blocks did not converge, and the length and range of each branch.
# The failed Student-t blocks are the missing dates of the common sample.
print_rolling_summary <- function(rolling) {
  cat("\n--- Rolling engine ---\n")
  cat("blocks:", rolling$n_blocks, "\n")
  cat("failed, Gaussian branch:", length(rolling$failed_normal),
      "| blocks:", rolling$failed_normal, "\n")
  cat("failed, Student-t branch:", length(rolling$failed_t),
      "| blocks:", rolling$failed_t, "\n")
  cat("rows, Gaussian branch:", nrow(rolling$normal),
      "| range:", as.character(range(rolling$normal$date)), "\n")
  cat("rows, Student-t branch:", nrow(rolling$t),
      "| range:", as.character(range(rolling$t$date)), "\n")
}


# Split one branch into its CCC and DCC series, with the same columns
# (date, sigma_p, sigma_spy, sigma_tlt): two cells of the 2x2 design.
split_ccc_dcc <- function(branch) {
  ccc <- data.frame(date = branch$date,
                    sigma_p = branch$sigma_p_ccc,
                    sigma_spy = branch$sigma_spy,
                    sigma_tlt = branch$sigma_tlt)
  dcc <- data.frame(date = branch$date,
                    sigma_p = branch$sigma_p_dcc,
                    sigma_spy = branch$sigma_spy,
                    sigma_tlt = branch$sigma_tlt)
  list(ccc = ccc, dcc = dcc)
}


# 6. From portfolio sigma to VaR

# Attach the realized portfolio return to each forecast date and standardize
# it by the forecast sigma.
add_standardized_returns <- function(sigmas, returns) {
  merged <- merge(sigmas, returns[, c("Date", "r_6040")],
                  by.x = "date", by.y = "Date", all.x = TRUE)
  merged$z_p <- merged$r_6040 / merged$sigma_p
  merged
}


# One-day VaR, reported as a positive loss: VaR = -(mu_p + q_alpha * sigma_p),
# with mu_p = 0 (paper, eq. 3). A violation occurs when r_p < -VaR.

# Gaussian tail: fixed quantile, -1.645 at 5% and -2.326 at 1%.
gaussian_var <- function(sigma_p, alpha) {
  -qnorm(alpha) * sigma_p
}


# Student-t tail at portfolio level: the degrees of freedom of a standardized
# t (unit variance) are estimated on the previous tail_window standardized
# returns and refreshed every refit_every days, as for the volatility models.
student_t_var <- function(z_p, sigma_p, alpha, tail_window, refit_every) {
  n <- length(z_p)
  var_values <- rep(NA_real_, n)
  nu <- NA_real_
  for (day in (tail_window + 1):n) {
    is_refit_day <- (day - tail_window - 1) %% refit_every == 0
    if (is_refit_day || is.na(nu)) {
      past_z <- z_p[(day - tail_window):(day - 1)]
      nu_estimate <- tryCatch(stdFit(past_z)$par["nu"], error = function(e) NA)
      if (!is.na(nu_estimate)) {
        nu <- nu_estimate   # if the fit fails, the previous value is kept
      }
    }
    if (is.na(nu)) next
    q_alpha <- qstd(alpha, mean = 0, sd = 1, nu = nu)
    var_values[day] <- -q_alpha * sigma_p[day]
  }
  var_values
}


# Semiparametric tail (filtered historical simulation; Barone-Adesi et al.,
# 1999): the empirical alpha-quantile of the previous tail_window standardized
# returns, rescaled by today's sigma. No distributional shape is imposed.
semiparametric_var <- function(z_p, sigma_p, alpha, tail_window) {
  n <- length(z_p)
  var_values <- rep(NA_real_, n)
  for (day in (tail_window + 1):n) {
    past_z <- z_p[(day - tail_window):(day - 1)]
    q_alpha <- quantile(past_z, probs = alpha, names = FALSE)
    var_values[day] <- -q_alpha * sigma_p[day]
  }
  var_values
}


# Historical-simulation benchmark: the empirical alpha-quantile of the previous
# hs_window raw portfolio returns. No volatility model at all.
historical_var <- function(all_returns, all_dates, target_dates, alpha, hs_window) {
  positions <- match(target_dates, all_dates)
  var_values <- rep(NA_real_, length(target_dates))
  for (i in seq_along(positions)) {
    past_returns <- all_returns[(positions[i] - hs_window):(positions[i] - 1)]
    var_values[i] <- -quantile(past_returns, probs = alpha, names = FALSE)
  }
  var_values
}


# VaR table for one level alpha, one row per date of the Gaussian branch. The
# Student-t branch misses the dates of its non-converged blocks, so its
# columns are joined by date and left NA there.
var_table <- function(alpha, normal_ccc, normal_dcc, t_ccc, t_dcc,
                      returns, tail_window, refit_every, hs_window) {
  var_values <- data.frame(date = normal_dcc$date, r_6040 = normal_dcc$r_6040)
  var_values$CCC_N <- gaussian_var(normal_ccc$sigma_p, alpha)
  var_values$DCC_N <- gaussian_var(normal_dcc$sigma_p, alpha)

  t_columns <- data.frame(
    date = t_ccc$date,
    CCC_t = student_t_var(t_ccc$z_p, t_ccc$sigma_p, alpha, tail_window, refit_every),
    DCC_t = student_t_var(t_dcc$z_p, t_dcc$sigma_p, alpha, tail_window, refit_every),
    CCC_fhs = semiparametric_var(t_ccc$z_p, t_ccc$sigma_p, alpha, tail_window),
    DCC_fhs = semiparametric_var(t_dcc$z_p, t_dcc$sigma_p, alpha, tail_window)
  )
  var_values <- merge(var_values, t_columns, by = "date", all.x = TRUE)

  var_values$HS <- historical_var(returns$r_6040, returns$Date, var_values$date, alpha, hs_window)
  var_values
}


# 7. Figure 3 

# Forecast 5% VaR against realized daily losses, October 2021 to March 2023.
plot_var_2022 <- function(var_values, path) {
  models <- c("CCC_N", "DCC_N", "CCC_t", "DCC_t", "HS")
  colors <- c(CCC_N = "blue", DCC_N = "orange", CCC_t = "green4",
              DCC_t = "red", HS = "grey40")

  in_window <- var_values$date >= as.Date("2021-10-01") & var_values$date <= as.Date("2023-03-31")
  zoom <- var_values[in_window, ]
  zoom$loss <- -zoom$r_6040

  png(path, width = 1200, height = 500)
  plot(zoom$date, zoom$loss, type = "h", col = "grey80",
       xlab = "", ylab = "VaR / perdita (giornaliero)",
       main = "VaR 5% vs perdita - zoom 2022",
       ylim = range(c(zoom$loss, zoom[models]), na.rm = TRUE))
  abline(h = 0, col = "black", lwd = 0.5)
  for (model in models) {
    lines(zoom$date, zoom[[model]], col = colors[model], lwd = 1.3)
  }
  legend("topleft", legend = c("Perdita realizzata", models),
         col = c("grey80", colors[models]), lwd = 2, ncol = 3, cex = 0.9)
  invisible(dev.off()) #close and save without printing
}


# Main
dir.create(RESULTS_DIR, showWarnings = FALSE)
dir.create(FIGURES_DIR, showWarnings = FALSE)

# 1) Data and specifications
inputs <- load_returns(INPUT_FILE)
spec_normal <- marginal_spec("norm")
spec_t <- marginal_spec("std")
dcc_normal <- dcc_spec(spec_normal)
dcc_t <- dcc_spec(spec_t)

# 2) Full-sample estimation: diagnostics on screen, Table 6 to file
marginal_fits <- fit_marginals_full_sample(inputs$asset_returns, spec_normal, spec_t)
print_marginal_diagnostics(marginal_fits)
marginal_parameters <- marginal_parameter_table(marginal_fits)
dcc_parameters <- dcc_parameter_table(dcc_normal, dcc_t, inputs$asset_returns)
print(marginal_parameters)
print(dcc_parameters)
write.csv(marginal_parameters, file.path(RESULTS_DIR, "params_marginals.csv"))
write.csv(dcc_parameters, file.path(RESULTS_DIR, "params_dcc.csv"))

# 3) Rolling engine: portfolio sigma forecasts for both branches
rolling <- rolling_engine(dcc_normal, dcc_t, inputs$asset_returns, WEIGHTS,
                          WINDOW_SIZE, REFIT_EVERY)
print_rolling_summary(rolling)
normal_branch <- split_ccc_dcc(rolling$normal)
t_branch <- split_ccc_dcc(rolling$t)

write.csv(normal_branch$ccc, file.path(RESULTS_DIR, "sigma_ccc_norm.csv"), row.names = FALSE)
write.csv(normal_branch$dcc, file.path(RESULTS_DIR, "sigma_dcc_norm.csv"), row.names = FALSE)
write.csv(t_branch$ccc, file.path(RESULTS_DIR, "sigma_ccc_t.csv"), row.names = FALSE)
write.csv(t_branch$dcc, file.path(RESULTS_DIR, "sigma_dcc_t.csv"), row.names = FALSE)

# 4) VaR at 5% and 1%
normal_ccc <- add_standardized_returns(normal_branch$ccc, inputs$returns)
normal_dcc <- add_standardized_returns(normal_branch$dcc, inputs$returns)
t_ccc <- add_standardized_returns(t_branch$ccc, inputs$returns)
t_dcc <- add_standardized_returns(t_branch$dcc, inputs$returns)

var_5pct <- var_table(0.05, normal_ccc, normal_dcc, t_ccc, t_dcc, inputs$returns,
                      TAIL_WINDOW, REFIT_EVERY, HS_WINDOW)
var_1pct <- var_table(0.01, normal_ccc, normal_dcc, t_ccc, t_dcc, inputs$returns,
                      TAIL_WINDOW, REFIT_EVERY, HS_WINDOW)

write.csv(var_5pct, file.path(RESULTS_DIR, "var_5pct.csv"), row.names = FALSE)
write.csv(var_1pct, file.path(RESULTS_DIR, "var_1pct.csv"), row.names = FALSE)

# 5) Figure 3
plot_var_2022(var_5pct, file.path(FIGURES_DIR, "rep_var5_2022.png"))

cat("[OK] tables in", RESULTS_DIR, "| figure in", FIGURES_DIR, "\n")
