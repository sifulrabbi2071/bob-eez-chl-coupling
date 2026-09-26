"""
============================================================================
04_coupling_analysis.py
----------------------------------------------------------------------------
Purpose : Seasonal partial correlations between deseasonalized surface Chl-a
          anomalies and each physical driver (SST, SSS), controlling for the
          other; plus two robustness tests: (i) months with Chl-a coverage
          below 15% excluded, with 10% and 20% as sensitivity thresholds
          (Table S2); (ii) the same partial correlations on the
          21 seasonal means, one value per year, which removes serial
          dependence between months (Table S4, part A).
Input   : BoB_EEZ_master_with_anomalies.csv   (from script 02)
Output  : TableS2_seasonal_coupling.csv (full record and 10/15/20% thresholds),
          TableS3_chl_coverage.csv,
          TableS4a_coupling_seasonal_means.csv  + console summary
Method  : Partial correlation is computed from the residuals of each variable
          regressed on the controlling variable; significance from a two-tailed
          t-test on n-3 degrees of freedom.
============================================================================
"""
import numpy as np
import pandas as pd
from scipy import stats

M_IN    = "BoB_EEZ_master_with_anomalies.csv"
CSV_OUT = "TableS2_seasonal_coupling.csv"

m = pd.read_csv(M_IN)
SEASONS = ["pre_monsoon", "sw_monsoon", "post_monsoon", "winter"]

def star(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"

def partial_corr(x, y, z):
    """Partial correlation of x and y controlling for z, with two-tailed p."""
    def resid(a, b):
        b = np.column_stack([np.ones_like(b), b])
        beta = np.linalg.lstsq(b, a, rcond=None)[0]
        return a - b @ beta
    r  = np.corrcoef(resid(x, z), resid(y, z))[0, 1]
    dof = len(x) - 3
    t  = r * np.sqrt(dof / (1 - r**2))
    p  = 2 * stats.t.sf(abs(t), dof)
    return r, p

# full dataset and copies with low-coverage months excluded
THRESHOLDS = [10, 15, 20]          # 15% = main analysis; 10% and 20% = sensitivity
COLS = ["sst_C_anom", "sss_PSU_anom", "chl_log_anom"]
masked = {}
for th in THRESHOLDS:
    mt = m.copy()
    mt.loc[mt["chl_validpct"] < th, ["chl_geomean", "chl_log_anom"]] = np.nan
    masked[th] = mt

print("=" * 78)
print("SEASONAL PARTIAL CORRELATION with surface Chl-a  (full | thresholds 10/15/20%)")
print("=" * 78)
rows = []
for s in SEASONS:
    d0 = m[m["season"] == s].dropna(subset=COLS)
    for a, b, ctrl, name in [
            ("sst_C_anom",   "chl_log_anom", "sss_PSU_anom", "SST-Chl (control SSS)"),
            ("sss_PSU_anom", "chl_log_anom", "sst_C_anom",   "SSS-Chl (control SST)")]:
        r0, p0 = partial_corr(d0[a].values, d0[b].values, d0[ctrl].values)
        row = {"Season": s, "Coupling": name,
               "partial_r_full": round(r0, 2), "p_full": round(p0, 3),
               "sig_full": star(p0), "n_full": len(d0)}
        txt = f"  {s:13s} {name:22s}: full {r0:+.2f}({star(p0)})"
        for th in THRESHOLDS:
            mt = masked[th]
            d1 = mt[mt["season"] == s].dropna(subset=COLS)
            r1, p1 = partial_corr(d1[a].values, d1[b].values, d1[ctrl].values)
            row.update({f"partial_r_thresh{th}": round(r1, 2), f"p_thresh{th}": round(p1, 3),
                        f"sig_thresh{th}": star(p1), f"n_thresh{th}": len(d1)})
            txt += f" | {th}% {r1:+.2f}({star(p1)})"
        rows.append(row)
        print(txt)

pd.DataFrame(rows).to_csv(CSV_OUT, index=False)
print(f"\nSaved {CSV_OUT}")

# ---- seasonal-mean test (one value per year; Table S4, part A) ----
print("\n" + "=" * 78)
print("SEASONAL-MEAN TEST  (partial correlation on 21 yearly seasonal means)")
print("=" * 78)
rows = []
for s in SEASONS:
    y = (m[m["season"] == s].groupby("year")[["sst_C_anom", "sss_PSU_anom", "chl_log_anom"]]
         .mean().dropna())
    for a, ctrl, name in [("sst_C_anom", "sss_PSU_anom", "SST-Chl (control SSS)"),
                          ("sss_PSU_anom", "sst_C_anom", "SSS-Chl (control SST)")]:
        r, p = partial_corr(y[a].values, y["chl_log_anom"].values, y[ctrl].values)
        rows.append({"Season": s, "Coupling": name, "partial_r": round(r, 2),
                     "p": round(p, 3), "sig": star(p), "n_years": len(y)})
        print(f"  {s:13s} {name:22s}: {r:+.2f} ({star(p)}, p={p:.3f}, n={len(y)})")
pd.DataFrame(rows).to_csv("TableS4a_coupling_seasonal_means.csv", index=False)
print("Saved TableS4a_coupling_seasonal_means.csv")

# coverage summary (Table S3)
cov = m.groupby("season")["chl_validpct"].agg(["mean", "min", "max"]).round(1).reindex(SEASONS)
cov["months_below_15pct"] = [int((m[m["season"] == s]["chl_validpct"] < 15).sum()) for s in SEASONS]
cov.to_csv("TableS3_chl_coverage.csv")
print("Saved TableS3_chl_coverage.csv")
print("\nChl-a valid-pixel coverage (%):")
print(cov.to_string())
