# Functions shared by the VaR backtest (03_backtest.R) and the robustness
# checks (04_robustness.R). This file only defines functions and constants;

# Tests used:
#   - unconditional coverage (Kupiec 1995)
#   - independence and conditional coverage (Christoffersen 1998)
#   - Basel traffic light (BCBS 1996)

library(rugarch)   # VaRTest: unconditional and conditional coverage test

STUDY_PERIODS <- list(
  full_oos      = c("2009-12-16", "2024-12-30"), # the whole common out-of-sample period.
  year_2022     = c("2022-01-01", "2022-12-31"),# the 2022 calendar year
  post_boundary = c("2022-11-10", "2024-12-30") # the new correlation regime, once it has settled.
)



# Load the VaR series and keep only the common sample
load_common_sample <- function(var_file) {
  var_data <- read.csv(var_file, stringsAsFactors = FALSE)
  var_data$date <- as.Date(var_data$date)
  core_models <- c("CCC_N", "DCC_N", "CCC_t", "DCC_t")
  all_core_available <- complete.cases(var_data[, core_models])
  common_sample <- var_data[all_core_available, ]
  return(common_sample)
}


# Keep the rows between two dates (both included).
subset_period <- function(var_data, start_date, end_date) {
  in_period <- var_data$date >= start_date & var_data$date <= end_date
  period_data <- var_data[in_period, ]
  return(period_data)
}



# Violation sequence
violation_sequence <- function(returns, var_loss) {
  violations <- as.integer(returns < -var_loss)
  return(violations)
}



# Coverage tests computed directly from the violation sequence
# Used for the independence statistic (VaRTest does not report it separately)
# and as a fallback where VaRTest fails or returns NaN.
bernoulli_loglik <- function(n_zeros, n_ones, prob) {
  loglik_zeros <- 0
  if (n_zeros > 0) {
    loglik_zeros <- n_zeros * log(1 - prob)
  }
  loglik_ones <- 0
  if (n_ones > 0) {
    loglik_ones <- n_ones * log(prob)
  }
  return(loglik_zeros + loglik_ones)
}


christoffersen_tests <- function(violations, p) {
  n_days <- length(violations)
  n_violations <- sum(violations)

  # unconditional coverage 
  observed_rate <- n_violations / n_days
  loglik_nominal <- bernoulli_loglik(n_days - n_violations, n_violations, p)
  loglik_observed <- bernoulli_loglik(n_days - n_violations, n_violations, observed_rate)
  uc_LR <- -2 * (loglik_nominal - loglik_observed)
  uc_pval <- 1 - pchisq(uc_LR, df = 1)

  #count transitions between consecutive days
  yesterday <- violations[-n_days]
  today <- violations[-1]
  n00 <- sum(yesterday == 0 & today == 0)
  n01 <- sum(yesterday == 0 & today == 1)
  n10 <- sum(yesterday == 1 & today == 0)
  n11 <- sum(yesterday == 1 & today == 1)

  # The independence test compares what happens after a violation with what
  # happens after a quiet day.
  no_day_after_violation <- (n10 + n11) == 0
  no_day_after_quiet <- (n00 + n01) == 0
  if (no_day_after_violation || no_day_after_quiet) {
    result <- list(uc_LR = uc_LR, uc_pval = uc_pval,
                   ind_LR = NA, ind_pval = NA,
                   cc_LR = NA, cc_pval = NA)
    return(result)
  }

  # independence (Christoffersen 1998) 
  prob_any_day <- (n01 + n11) / (n00 + n01 + n10 + n11)
  prob_after_quiet <- n01 / (n00 + n01)
  prob_after_violation <- n11 / (n10 + n11)
  loglik_one_prob <- bernoulli_loglik(n00 + n10, n01 + n11, prob_any_day)
  loglik_after_quiet <- bernoulli_loglik(n00, n01, prob_after_quiet)
  loglik_after_violation <- bernoulli_loglik(n10, n11, prob_after_violation)
  loglik_two_probs <- loglik_after_quiet + loglik_after_violation
  ind_LR <- -2 * (loglik_one_prob - loglik_two_probs)
  ind_pval <- 1 - pchisq(ind_LR, df = 1)

  # conditional coverage 
  cc_LR <- uc_LR + ind_LR
  cc_pval <- 1 - pchisq(cc_LR, df = 2)

  result <- list(uc_LR = uc_LR, uc_pval = uc_pval,
                 ind_LR = ind_LR, ind_pval = ind_pval,
                 cc_LR = cc_LR, cc_pval = cc_pval)
  return(result)
}



# Coverage tests: rugarch first, above computation where it fails
# VaRTest (rugarch) gives unconditional and conditional coverage.
# The independence statistic always comes from our own computation, because
# VaRTest does not report it separately.
coverage_tests <- function(returns, var_loss, p) {
  violations <- violation_sequence(returns, var_loss)
  own_tests <- christoffersen_tests(violations, p)

  # rugarch expects the VaR as a return quantile (a negative number),
  # not as a positive loss.
  var_quantile <- -var_loss
  rugarch_result <- tryCatch(
    VaRTest(alpha = p, actual = returns, VaR = var_quantile, conf.level = 0.95),
    error = function(e) NULL
  )

  rugarch_failed <- is.null(rugarch_result)
  if (!rugarch_failed) {
    rugarch_failed <- !is.finite(rugarch_result$uc.LRstat) || !is.finite(rugarch_result$cc.LRstat)
  }
  if (rugarch_failed) {
    return(own_tests)
  }

  result <- list(uc_LR = rugarch_result$uc.LRstat,
                 uc_pval = rugarch_result$uc.LRp,
                 ind_LR = own_tests$ind_LR,
                 ind_pval = own_tests$ind_pval,
                 cc_LR = rugarch_result$cc.LRstat,
                 cc_pval = rugarch_result$cc.LRp)
  return(result)
}



# Basel traffic-light zone (BCBS 1996)

# The Basel rule is written for 250 days at the 1% level: green with 0-4
# violations, yellow with 5-9, red with 10 or more. 
basel_zone <- function(n_violations, n_days, p) {
  prob_this_many_or_fewer <- pbinom(n_violations, size = n_days, prob = p)
  if (prob_this_many_or_fewer <= 0.95) {
    return("green")
  }
  if (prob_this_many_or_fewer <= 0.9999) {
    return("yellow")
  }
  return("red")
}



# Backtest of one VaR series: one row of a results table
# Rounding as in the paper tables: expected 2 decimals, rate 4, LR 3, p-value 4.
backtest_one_series <- function(period_name, model_name, level_label, returns, var_loss, p) {
  violations <- violation_sequence(returns, var_loss)
  n_days <- length(violations)
  n_violations <- sum(violations)
  tests <- coverage_tests(returns, var_loss, p)

  result_row <- data.frame(
    period     = period_name,
    model      = model_name,
    level      = level_label,
    n_obs      = n_days,
    expected   = round(p * n_days, 2),
    actual     = n_violations,
    rate       = round(n_violations / n_days, 4),
    uc_LR      = round(tests$uc_LR, 3),
    uc_pval    = round(tests$uc_pval, 4),
    ind_LR     = round(tests$ind_LR, 3),
    ind_pval   = round(tests$ind_pval, 4),
    cc_LR      = round(tests$cc_LR, 3),
    cc_pval    = round(tests$cc_pval, 4),
    basel_zone = basel_zone(n_violations, n_days, p),
    stringsAsFactors = FALSE
  )
  return(result_row)
}



# Results table: every model, both levels, every period

# Row order: for each period, the models at 1%, then the models at 5%.
backtest_table_by_period <- function(var_1pct, var_5pct, periods, models) {
  rows <- list()
  for (period_name in names(periods)) {
    start_date <- as.Date(periods[[period_name]][1])
    end_date <- as.Date(periods[[period_name]][2])
    data_1pct <- subset_period(var_1pct, start_date, end_date)
    data_5pct <- subset_period(var_5pct, start_date, end_date)

    for (model in models) {
      row <- backtest_one_series(period_name, model, "1%",
                                 data_1pct$r_6040, data_1pct[[model]], 0.01)
      rows[[length(rows) + 1]] <- row
    }
    for (model in models) {
      row <- backtest_one_series(period_name, model, "5%",
                                 data_5pct$r_6040, data_5pct[[model]], 0.05)
      rows[[length(rows) + 1]] <- row
    }
  }
  # stack the one-row tables into a single table
  results_table <- do.call(rbind, rows)
  return(results_table)
}
