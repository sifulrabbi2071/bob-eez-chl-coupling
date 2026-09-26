"""
============================================================================
03b_seasonal_autocorrelation.py
----------------------------------------------------------------------------
Purpose : Check for inter-annual serial dependence in the per-season trend
          tests of script 03 (Mann-Kendall on the 21 yearly seasonal means).
          For each variable and season, the Sen's-slope trend is removed and
          the lag-1 autocorrelation of the residual series is compared with
          the approximate 95% limit (+/- 1.96 / sqrt(n)).
Input   : BoB_EEZ_master_with_anomalies.csv   (from script 02)
Output  : TableS1b_seasonal_autocorrelation.csv  + console summary
============================================================================
"""
import numpy as np
import pandas as pd
import pymannkendall as mk

m = pd.read_csv("BoB_EEZ_master_with_anomalies.csv")
VARS = [("sst_C", "SST"), ("sss_PSU", "SSS"), ("chl_geomean", "Chl-a")]
SEASONS = ["pre_monsoon", "sw_monsoon", "post_monsoon", "winter"]

rows = []
for col, lab in VARS:
    for s in SEASONS:
        y = m[m["season"] == s].dropna(subset=[col]).groupby("year")[col].mean().values
        r = mk.original_test(y)
        resid = y - (r.intercept + r.slope * np.arange(len(y)))   # remove Sen's-slope trend
        ac1 = np.corrcoef(resid[:-1], resid[1:])[0, 1]
        lim = 1.96 / np.sqrt(len(y))
        rows.append({"Variable": lab, "Season": s, "n_years": len(y),
                     "lag1_autocorr_detrended": round(ac1, 2), "approx_95pct_limit": round(lim, 2),
                     "significant": abs(ac1) > lim})
        print(f"  {lab:6s} {s:13s}: r1 = {ac1:+.2f}  (limit +/-{lim:.2f})")
out = pd.DataFrame(rows)
out.to_csv("TableS1b_seasonal_autocorrelation.csv", index=False)
print(f"\nmax |r1| = {out['lag1_autocorr_detrended'].abs().max():.2f}; "
      f"any significant: {bool(out['significant'].any())}")
print("Saved TableS1b_seasonal_autocorrelation.csv")
