# 60/40 portfolio VaR across the 2022 stock-bond regime shift

Code and results for the paper *The Impact of the 2022 Stock–Bond Regime Shift on 60/40 Daily Value-at-Risk* (Francesco Steffenino, 2026), written for the financial econometrics course project. The paper is in `paper/`.

## Question

In 2022 the correlation between US stocks and long-term Treasuries turned from negative to positive, and a 60/40 portfolio lost most of its diversification. The paper asks whether a dynamic correlation (DCC instead of CCC) or fat-tailed innovations (Student-t instead of Gaussian) improve the calibration of one-day VaR for a 60/40 portfolio across that shift, and which of the two matters more.

## Data

Daily adjusted closes of SPY (S&P 500) and TLT (20+ year Treasuries) from Yahoo Finance, January 2004 to December 2024: 5,284 daily log returns. The portfolio return is 0.6 · r_SPY + 0.4 · r_TLT with fixed weights. The regime-shift boundary, the first date after which the 252-day rolling correlation stays positive, is 10 November 2022.

## Method

- Marginals: AR(1)-GARCH(1,1) for each asset, with Gaussian or Student-t innovations.
- Correlation: CCC or DCC(1,1) on the standardized residuals. The 2 × 2 grid isolates the effect of each feature.
- Estimation: rolling 1,000-day window, re-estimated every 21 days, one-step-ahead forecasts of the portfolio volatility.
- VaR at 5% and 1% with Gaussian, Student-t and semiparametric (empirical quantile of the filtered residuals) tails, plus rolling historical simulation (250 days) as benchmark.
- Backtest on a common sample of 3,616 days (16 December 2009 to 30 December 2024): unconditional coverage (Kupiec), independence and conditional coverage (Christoffersen), Basel traffic light. Three periods: full sample, 2022, post-boundary.

## Main results

- Tail calibration matters most. Parametric tails are too conservative (2.8–3.8% violations at the 5% level); the semiparametric tail is close to nominal (5.2%).
- Volatility dynamics come next. Historical simulation has the right average rate (5.6%) but clustered violations, and in 2022 it falls in the Basel red zone (10.8%).
- Correlation is a second-order effect. In 2022 the DCC is near nominal (5.6%) while the CCC is too permissive (7.2–7.6%, yellow zone). The advantage does not last: after the boundary the DCC becomes conservative (2.7–2.9%), and with the semiparametric tail the 2022 gap shrinks to 6.0% against 6.4%.

At a one-day horizon, portfolio variance is dominated by the two volatilities, so even a large change in correlation moves the VaR little.

## Repository

```
src/
  00_download_data.py     prices from Yahoo Finance -> data/
  01_data_analysis.py     returns, descriptive statistics, rolling correlation,
                          variance decomposition, Ljung-Box (Section 3)
  02_estimate_var.R       GARCH/CCC/DCC, rolling engine, VaR (Section 4)
  backtest_functions.R    coverage tests and Basel zones, used by 03 and 04
  03_backtest.R           backtest by period (Section 5)
  04_robustness.R         semiparametric tail, boundary shifted by ±3 months
results/                  CSV tables
figures/                  figures of the paper
paper/                    the paper (PDF)
requirements.txt          Python packages
R_packages.txt            R session and package versions
```

| Paper | File |
|---|---|
| Table 1 | `results/descriptive_stats.csv` |
| Table 2 | `results/variance_decomposition.csv`, `results/drawdown_2022.csv` |
| Table 3 | `results/ljungbox_meaneq.csv` |
| Table 4 | `results/backtest_by_period.csv` |
| Table 5 | `results/backtest_fhs.csv` |
| Table 6 | `results/params_marginals.csv`, `results/params_dcc.csv` |
| Table 7 | `results/boundary_robustness.csv` |
| Figures 1, 2, 3 | `figures/correlation_63_252.png`, `figures/acf_60_40.png`, `figures/rep_var5_2022.png` |

## How to reproduce

Requirements: Python 3 with the packages in `requirements.txt`, R with rugarch, rmgarch, xts and fGarch. Run everything from the repository root:

```
python src/00_download_data.py
python src/01_data_analysis.py
Rscript src/02_estimate_var.R 
Rscript src/03_backtest.R
Rscript src/04_robustness.R
```

Each step reads the CSVs written by the previous one, so Step 3 and 4 can be run alone on the files already in `results/`.

## Reproducibility notes

- Prices. Yahoo revises adjusted closes over time, so a new download can differ slightly from the prices used in the paper. The files in `results/` were produced with the prices downloaded at the time of the analysis.
- Student-t branch. In a few estimation windows (2013, 2015, 2024; none in 2022) the default optimizer fails and a fallback with random restarts is used. A new run can change which of these windows converge and slightly alter the Student-t and semiparametric results in those years. The CSVs in `results/` are the ones behind the paper.
- Coverage tests. The conditional-coverage statistic for historical simulation at 5% on the full sample is computed directly from the violation sequence, because rugarch's `VaRTest` returns NaN there. For the same numerical reason (very small likelihoods), `VaRTest` can return the unconditional-coverage p-value of the two semiparametric 5% full-sample rows with a small difference in the fourth decimal depending on the platform. These p-values are not reported in the paper.
