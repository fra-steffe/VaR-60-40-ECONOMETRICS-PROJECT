#main

# Backtest of the one-day portfolio VaR: four core models (CCC/DCC x
# Gaussian/Student-t) and the Historical Simulation benchmark, at the 1% and
# 5% levels, over three periods (full sample, 2022, post-boundary).

source("src/backtest_functions.R")


var_1pct_file <- "results/var_1pct.csv"
var_5pct_file <- "results/var_5pct.csv"
output_file <- "results/backtest_by_period.csv"

# 1. Load both VaR files on the common sample (see load_common_sample).
var_1pct <- load_common_sample(var_1pct_file)
var_5pct <- load_common_sample(var_5pct_file)

# The two files must cover exactly the same days, or the 1% and 5% results
# would not be comparable.
stopifnot(identical(var_1pct$date, var_5pct$date))
cat("Common sample:", nrow(var_1pct), "days, from",
    format(min(var_1pct$date)), "to", format(max(var_1pct$date)), "\n")

# 2. Backtest every model, at both levels, in every period.
models <- c("CCC_N", "DCC_N", "CCC_t", "DCC_t", "HS")
results <- backtest_table_by_period(var_1pct, var_5pct, STUDY_PERIODS, models)

# 3. Save. Statistics that are not defined (no violations) are left empty.
write.csv(results, output_file, row.names = FALSE, na = "")
print(results, row.names = FALSE)

# 4. Check the cell that VaRTest cannot compute: HS at 5% on the full sample.
#    Expected: cc_LR = 22.805, cc_pval below 0.001.
is_hs_5pct_full <- results$period == "full_oos" & results$model == "HS" & results$level == "5%"
cat("\nHS, 5%, full sample (computed directly, VaRTest returns NaN here):\n")
print(results[is_hs_5pct_full, c("uc_LR", "ind_LR", "cc_LR", "cc_pval")], row.names = FALSE)
