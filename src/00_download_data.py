
#  Step 0 - Download daily prices of SPY and TLT from Yahoo Finance.
import os

import pandas as pd
import yfinance as yf

OUTPUT_FILE = "data/prices_spy_tlt.csv"
TICKERS = ["SPY", "TLT"]
START_DATE = "2004-01-01"
END_DATE = "2025-01-01"   # exclusive: the last day downloaded is 2024-12-31


def download_adjusted_closes():
    # auto_adjust=False keeps the "Adj Close" column separate from "Close".
    raw = yf.download(TICKERS, start=START_DATE, end=END_DATE,
                      auto_adjust=False, progress=False)
    adjusted = raw["Adj Close"][TICKERS]
    prices = pd.DataFrame({
        "SPY_AdjClose": adjusted["SPY"],
        "TLT_AdjClose": adjusted["TLT"],
    })
    prices.index.name = "Date"
    prices = prices.dropna()
    return prices


def main():
    os.makedirs("data", exist_ok=True)
    prices = download_adjusted_closes()
    prices.to_csv(OUTPUT_FILE)
    print("[data]", len(prices), "days,",
          prices.index.min().date(), "->", prices.index.max().date())
    print("[OK] saved", OUTPUT_FILE)


if __name__ == "__main__":
    main()
