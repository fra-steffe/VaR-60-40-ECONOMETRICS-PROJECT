
# Two robustness checks on the backtest results:
#   1. tail: the parametric Student-t quantile replaced by the semiparametric
#      one (empirical quantile of the standardised residuals, fhs)
#   2. boundary: the start of the post-boundary period moved by up to three
#      months around 2022-11-10.
source("src/backtest_functions.R")

# Tail robustness: Student-t quantile vs semiparametric quantile

# The two versions share the same GARCH-t marginals and the same correlation
# (CCC or DCC): only the tail quantile changes.
tail_robustness <- function(var_1pct, var_5pct) {
  tail_models <- c("CCC_t", "DCC_t", "CCC_fhs", "DCC_fhs")
  full_table <- backtest_table_by_period(var_1pct, var_5pct, STUDY_PERIODS, tail_models)
  columns_to_keep <- c("period", "model", "level", "n_obs", "expected", "actual",
                       "rate", "uc_pval", "cc_pval", "basel_zone")
  tail_table <- full_table[, columns_to_keep]
  return(tail_table)
}


# Boundary robustness: move the start of the post-boundary period

# we restart the period one, two and three months
# earlier and later, always ending on 2024-12-30, and repeat the 5% backtest.
# We use the Student-t arm, the one reported for this result in the paper.
boundary_robustness <- function(var_5pct, boundary_dates, end_date) {
  rows <- list()
  for (boundary in boundary_dates) {
    data_after_boundary <- subset_period(var_5pct, as.Date(boundary), as.Date(end_date))
    for (model in c("CCC_t", "DCC_t")) {
      backtest_row <- backtest_one_series(boundary, model, "5%",
                                          data_after_boundary$r_6040,
                                          data_after_boundary[[model]], 0.05)
      row <- data.frame(boundary = boundary,
                        model   = model,
                        n_obs   = backtest_row$n_obs,
                        actual  = backtest_row$actual,
                        rate    = backtest_row$rate,
                        cc_pval = backtest_row$cc_pval,
                        basel_zone = backtest_row$basel_zone,
                        stringsAsFactors = FALSE)
      rows[[length(rows) + 1]] <- row
    }
  }
  # stack the one-row tables into a single table
  boundary_table <- do.call(rbind, rows)
  return(boundary_table)
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

var_1pct_file <- "results/var_1pct.csv"
var_5pct_file <- "results/var_5pct.csv"
tail_output_file <- "results/backtest_fhs.csv"
boundary_output_file <- "results/boundary_robustness.csv"

# 1. Load both VaR files on the common sample.
var_1pct <- load_common_sample(var_1pct_file)
var_5pct <- load_common_sample(var_5pct_file)
stopifnot(identical(var_1pct$date, var_5pct$date))

# 2. Tail robustness.
tail_table <- tail_robustness(var_1pct, var_5pct)
write.csv(tail_table, tail_output_file, row.names = FALSE, na = "")
print(tail_table, row.names = FALSE)

# 3. Boundary robustness: seven start dates, one month apart, centred on
#    2022-11-10.
boundary_dates <- c("2022-08-10", "2022-09-10", "2022-10-10", "2022-11-10",
                    "2022-12-10", "2023-01-10", "2023-02-10")
boundary_table <- boundary_robustness(var_5pct, boundary_dates, "2024-12-30")
write.csv(boundary_table, boundary_output_file, row.names = FALSE, na = "")
print(boundary_table, row.names = FALSE)
