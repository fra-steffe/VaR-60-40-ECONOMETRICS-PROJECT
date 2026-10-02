#  Step 1 - Data and preliminary analysis (paper, Section 3).


import os

import numpy as np
import pandas as pd
import scipy
from scipy import stats

import statsmodels
from statsmodels.stats.diagnostic import acorr_ljungbox
from statsmodels.graphics.tsaplots import plot_acf

import matplotlib
matplotlib.use("Agg")  # backend senza finestra grafica: salva solo su file
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


# PARAMETERS
INPUT_FILE = "data/prices_spy_tlt.csv"   
RESULTS_DIR = "results"
FIGURES_DIR = "figures"
     

WEIGHT_SPY = 0.60 # fixed 60/40 weights, implicit daily rebalancing
WEIGHT_TLT = 0.40
SHORT_WINDOW = 63 # one quarter: timing of the correlation flip               
LONG_WINDOW = 252 # one year: level of the correlation regime               
LJUNG_BOX_LAGS = [5, 10, 20]        
TRADING_DAYS_PER_YEAR = 252         

# Load prices indexed by date, in chronological order
def load_prices(path):
    prices = pd.read_csv(path, parse_dates=["Date"])
    prices = prices.set_index("Date")
    prices = prices.sort_index()
    return prices


# Daily log returns on adjusted closes, so that TLT's coupon income is included.
# The portfolio return is the weighted sum of log returns (paper, eq 1)
def compute_log_returns(prices):
    spy_returns = np.log(prices["SPY_AdjClose"]).diff()
    tlt_returns = np.log(prices["TLT_AdjClose"]).diff()
    portfolio_returns = WEIGHT_SPY * spy_returns + WEIGHT_TLT * tlt_returns
    returns = pd.DataFrame({
        "r_SPY": spy_returns,
        "r_TLT": tlt_returns,
        "r_6040": portfolio_returns,
    })
    returns = returns.dropna() # first day has no previous price
    return returns


# Descriptive statistics of the three return series (paper, Table 1) 
def descriptive_stats(returns):
    stats_by_series = {}
    for column_name in returns.columns:
        series_values = returns[column_name].to_numpy()
        mean_daily = np.mean(series_values)
        sd_daily = np.std(series_values, ddof=1)
        stats_by_series[column_name] = {
            "n": series_values.size,
            "mean_daily": mean_daily,
            "mean_ann": mean_daily * TRADING_DAYS_PER_YEAR,
            "sd_daily": sd_daily,
            "sd_ann": sd_daily * np.sqrt(TRADING_DAYS_PER_YEAR), # volatility scales with sqrt(time)
            "min": np.min(series_values),
            "max": np.max(series_values),
            "skew": stats.skew(series_values),                     # skewness
            "excess_kurt": stats.kurtosis(series_values, fisher=True),  # excess kurtosis (gaussian = 0)
        }
    table = pd.DataFrame(stats_by_series).T
    return table


# Rolling SPY-TLT correlation on two windows: 63 daysand 252 days.
def rolling_correlation(returns):
    short_window_corr = returns["r_SPY"].rolling(SHORT_WINDOW).corr(returns["r_TLT"])
    long_window_corr = returns["r_SPY"].rolling(LONG_WINDOW).corr(returns["r_TLT"])
    correlations = pd.DataFrame({
        "corr63": short_window_corr, # first 62 rows are NaN
        "corr252": long_window_corr, # first 251 rows are NaN
    })
    return correlations


# Regime-shift boundary: the first date from which the 252-day correlation
# stays non-negative until the end of the sample (paper: 10 November 2022).
# Found as the day after the last negative value. Assumes the sample ends in
# the positive regime, which holds here.
def find_regime_boundary(rolling_correlations):
    long_corr = rolling_correlations["corr252"].dropna()
    is_negative = long_corr.to_numpy() < 0
    negative_positions = np.where(is_negative)[0]
    last_negative_position = negative_positions.max()
    boundary_date = long_corr.index[last_negative_position + 1]
    return boundary_date


# decomposes Portofolio variance in one specified sub-period, useful in 
# the function below. The portfolio variance is split into its three terms (paper, eq 2)
def decompose_one_period(returns, start, end, label):
    window = returns.loc[start:end]
    sd_spy = window["r_SPY"].std(ddof=1)
    sd_tlt = window["r_TLT"].std(ddof=1)
    correlation = window["r_SPY"].corr(window["r_TLT"])
    own_var_spy = (WEIGHT_SPY * sd_spy) ** 2
    own_var_tlt = (WEIGHT_TLT * sd_tlt) ** 2
    cross_term = 2 * WEIGHT_SPY * WEIGHT_TLT * correlation * sd_spy * sd_tlt
    total_var = own_var_spy + own_var_tlt + cross_term
    annualized_sd = np.sqrt(total_var) * np.sqrt(TRADING_DAYS_PER_YEAR)
    cross_share = cross_term / total_var * 100
    return {
        "period": label,
        "sd_SPY": sd_spy,
        "sd_TLT": sd_tlt,
        "rho": correlation,
        "var_own_SPY": own_var_spy,
        "var_own_TLT": own_var_tlt,
        "cross_term": cross_term,
        "var_total": total_var,
        "sd_port_ann": annualized_sd,
        "cross_pct_of_total": cross_share,
    }


# Repeat the decomposition over the four sub-periods of Table 2, to show the
# cross term changing sign between regimes.
def variance_decomposition(returns, regime_boundary):
    boundary_label = str(regime_boundary.date())
    periods = [
        ("2004", "2021-12-31", "old regime (2004-2021)"),
        ("2022-01-01", "2022-12-31", "2022"),
        (boundary_label, "2024-12-31", "regime shift (post-boundary)"),
        ("2004", "2024-12-31", "full sample"),
    ]
    decomposition_rows = []
    for start, end, label in periods:
        row = decompose_one_period(returns, start, end, label)
        decomposition_rows.append(row)
    table = pd.DataFrame(decomposition_rows)
    table = table.set_index("period")
    return table


# Maximum 2022 drawdown of the 60/40 portfolio
def compute_drawdown_2022(returns):
    portfolio_returns = returns["r_6040"]
    wealth_index = np.exp(portfolio_returns.cumsum())
    running_peak = wealth_index.cummax()
    drawdown = wealth_index / running_peak - 1

    drawdown_in_2022 = drawdown.loc["2022"]
    trough_date = drawdown_in_2022.idxmin()
    max_drawdown = drawdown_in_2022.min()
    peak_date = wealth_index.loc[:trough_date].idxmax()

    calendar_2022_return = np.exp(returns.loc["2022", "r_6040"].sum()) - 1

    summary = {
        "peak_date": peak_date.date(),
        "trough_date": trough_date.date(),
        "max_drawdown_2022": max_drawdown,
        "calendar_2022_return": calendar_2022_return,
    }
    return summary


# Ljung-Box tests for the mean equation (paper, Table 3). On return levels they
# decide between a constant mean and an AR(1); on squared returns they detect
# volatility clustering (ARCH effects), which motivates the GARCH models.
def mean_equation_diagnostics(returns):
    diagnostics_by_series = {}
    for column_name in returns.columns:
        series_values = returns[column_name].to_numpy()
        squared_values = series_values ** 2
        levels_test = acorr_ljungbox(series_values, lags=LJUNG_BOX_LAGS, return_df=True)
        squares_test = acorr_ljungbox(squared_values, lags=LJUNG_BOX_LAGS, return_df=True)

        row = {}
        row["acf1_lvl"] = returns[column_name].autocorr(lag=1)
        for lag in LJUNG_BOX_LAGS:
            row[f"LB{lag}_lvl"] = levels_test.loc[lag, "lb_stat"]
            row[f"LB{lag}_lvl_p"] = levels_test.loc[lag, "lb_pvalue"]
            row[f"LB{lag}_sq"] = squares_test.loc[lag, "lb_stat"]
            row[f"LB{lag}_sq_p"] = squares_test.loc[lag, "lb_pvalue"]
        diagnostics_by_series[column_name] = row

    table = pd.DataFrame(diagnostics_by_series).T
    return table


# Figure 1: 63-day (flip timing) and 252-day (regime level) correlation, with
# the regime-shift boundary.
def plot_correlation(rolling_correlations, regime_boundary, output_dir):
    short_corr = rolling_correlations["corr63"]
    long_corr = rolling_correlations["corr252"]
    boundary_label = "regime-shift boundary ({})".format(regime_boundary.date())

    figure, axis = plt.subplots(figsize=(9, 4.2))
    axis.axhline(0, color="0.5", lw=0.8)
    axis.plot(short_corr.index, short_corr.to_numpy(), color="#9ecae1", lw=1.3,
              label="63-day (flip timing)")
    axis.plot(long_corr.index, long_corr.to_numpy(), color="#1f4e79", lw=1.3,
              label="252-day (regime level)")
    axis.axvline(regime_boundary, color="#b22222", ls="--", lw=1.0, label=boundary_label)
    axis.set_ylabel("SPY-TLT correlation")
    axis.legend(loc="lower left", fontsize=8, frameon=False)
    axis.xaxis.set_major_locator(mdates.YearLocator(2))
    axis.margins(x=0.01)
    axis.set_title("SPY-TLT rolling correlation of daily returns (2004-2024)")
    figure.tight_layout()
    figure.savefig(output_dir + "/correlation_63_252.png", dpi=150)
    plt.close(figure)


# Figure 2: autocorrelation of 60/40 returns (left, almost flat) and squared
# returns (right, large and persistent: volatility clustering).
def plot_acf_charts(returns, output_dir):
    portfolio_returns = returns["r_6040"].to_numpy()
    squared_returns = portfolio_returns ** 2
    figure, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    plot_acf(portfolio_returns, lags=20, ax=axes[0], zero=False,
             title="ACF of 60/40 returns")
    plot_acf(squared_returns, lags=20, ax=axes[1], zero=False,
             title="ACF of 60/40 squared returns")
    figure.tight_layout()
    figure.savefig(output_dir + "/acf_60_40.png", dpi=150)
    plt.close(figure)


def main():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    os.makedirs(FIGURES_DIR, exist_ok=True)
 
    # 1) Prices and returns
    prices = load_prices(INPUT_FILE)
    returns = compute_log_returns(prices)
 
    # 2) Analysis
    descriptives = descriptive_stats(returns)
    rolling_correlations = rolling_correlation(returns)
    regime_boundary = find_regime_boundary(rolling_correlations)
    var_decomposition = variance_decomposition(returns, regime_boundary)
    drawdown_summary = compute_drawdown_2022(returns)
    ljung_box_table = mean_equation_diagnostics(returns)
 
    # 3) Quick check of the key numbers on screen (details are in the CSVs)
    print("[data]", len(returns), "daily returns,",
          returns.index.min().date(), "->", returns.index.max().date())
    print("[correlation] regime-shift boundary:", regime_boundary.date())
    print("[drawdown 2022] {:.4f} (peak {} -> trough {})".format(
        drawdown_summary["max_drawdown_2022"],
        drawdown_summary["peak_date"],
        drawdown_summary["trough_date"]))
 
    # 4) Tables
    returns.to_csv(RESULTS_DIR + "/returns_spy_tlt_6040.csv", float_format="%.10f")
    descriptives.to_csv(RESULTS_DIR + "/descriptive_stats.csv", float_format="%.6f")
    rolling_correlations.dropna().to_csv(RESULTS_DIR + "/rolling_correlation.csv", float_format="%.6f")
    var_decomposition.to_csv(RESULTS_DIR + "/variance_decomposition.csv", float_format="%.8f")
    pd.DataFrame([drawdown_summary]).to_csv(RESULTS_DIR + "/drawdown_2022.csv", index=False, float_format="%.6f")
    ljung_box_table.to_csv(RESULTS_DIR + "/ljungbox_meaneq.csv", float_format="%.4f")
 
    # 5) Figures
    plot_correlation(rolling_correlations, regime_boundary, FIGURES_DIR)
    plot_acf_charts(returns, FIGURES_DIR)
 
    print("[OK] tables in", RESULTS_DIR, "| figures in", FIGURES_DIR)
 
 
if __name__ == "__main__":
    main()
