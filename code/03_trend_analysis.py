"""
============================================================================
03_trend_analysis.py
----------------------------------------------------------------------------
Purpose : Long-term trends in SST, SSS and surface Chl-a.
          A. Full record: seasonal Mann-Kendall test in its Hirsch-Slack form
             (serial dependence between seasons; Hirsch & Slack 1984) with the
             seasonal Sen's slope  -> Table 2 of the paper.
             The uncorrected seasonal test is printed for transparency.
          B. Per season: Mann-Kendall + Sen's slope on the 21 annual values of
             each seasonal mean  -> Fig. 4 / Table S1.
          C. Sensitivity: continuous December-February winter.
Input   : BoB_EEZ_master_with_anomalies.csv   (from script 02)
Output  : Table2_annual_trends.csv, TableS1_seasonal_trends.csv,
          TableS1_winter_DecFeb_sensitivity.csv
Note    : Missing months are kept as NaN. Dropping them would shift every
          later month into the wrong calendar position in the seasonal test.
============================================================================
"""
import numpy as np
import pandas as pd
import pymannkendall as mk

M_IN = "BoB_EEZ_master_with_anomalies.csv"

m = pd.read_csv(M_IN).sort_values("date").reset_index(drop=True)
assert len(m) % 12 == 0, "the monthly record must contain whole years"
VARS = [("sst_C", "SST", "°C/decade"),
        ("sss_PSU", "SSS", "PSU/decade"),
        ("chl_geomean", "Chl-a", "mg m-3/decade")]
SEASONS = ["pre_monsoon", "sw_monsoon", "post_monsoon", "winter"]

def star(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"

# ---- A. full-record trend (Table 2) ----
print("=" * 78)
print("FULL-RECORD TREND  (seasonal Mann-Kendall, period = 12; NaN months kept)")
print("=" * 78)
rows = []
for col, lab, unit in VARS:
    x  = m[col].values                                  # keep NaN -> calendar alignment preserved
    hs = mk.correlated_seasonal_test(x, period=12)      # Hirsch-Slack (reported in the paper)
    pl = mk.seasonal_test(x, period=12)                 # uncorrected, for transparency
    rows.append({"Variable": lab, "Unit": unit, "Sen_slope_decade": round(hs.slope * 10, 4),
                 "p_HirschSlack": round(hs.p, 4), "Significance": star(hs.p),
                 "p_uncorrected": round(pl.p, 4), "n_valid_months": int(np.isfinite(x).sum())})
    print(f"  {lab:6s}: {hs.slope*10:+.3f} {unit}   Hirsch-Slack p={hs.p:.4f} ({star(hs.p)})"
          f"   [uncorrected p={pl.p:.4f}]")
pd.DataFrame(rows).to_csv("Table2_annual_trends.csv", index=False)

# ---- B. per-season trend (Fig. 4 / Table S1) ----
print("\n" + "=" * 78)
print("SEASONAL TREND  (Mann-Kendall + Sen's slope on yearly seasonal means)")
print("=" * 78)
rows = []
for col, lab, unit in VARS:
    for s in SEASONS:
        yr = m[m["season"] == s].dropna(subset=[col]).groupby("year")[col].mean()
        r  = mk.original_test(yr.values)
        rows.append({"Variable": lab, "Unit": unit, "Season": s,
                     "Sen_slope_decade": round(r.slope * 10, 4), "p_value": round(r.p, 4),
                     "Significance": star(r.p), "Trend": r.trend, "n_years": len(yr)})
        print(f"  {lab:6s} {s:13s}: {r.slope*10:+.4f} ({star(r.p)}, p={r.p:.3f})")
pd.DataFrame(rows).to_csv("TableS1_seasonal_trends.csv", index=False)

# ---- C. sensitivity: continuous Dec-Feb winter (December joins the following Jan-Feb) ----
print("\n" + "=" * 78)
print("SENSITIVITY  winter as a continuous Dec-Feb season (2006-2025 winters)")
print("=" * 78)
w = m[m["season"] == "winter"].copy()
w["winter_year"] = np.where(w["month"] == 12, w["year"] + 1, w["year"])
w = w[(w["winter_year"] >= 2006) & (w["winter_year"] <= 2025)]
rows = []
for col, lab, unit in VARS:
    yr = w.dropna(subset=[col]).groupby("winter_year")[col].mean()
    r  = mk.original_test(yr.values)
    rows.append({"Variable": lab, "Unit": unit, "Season": "winter_DecFeb_continuous",
                 "Sen_slope_decade": round(r.slope * 10, 4), "p_value": round(r.p, 4),
                 "Significance": star(r.p), "Trend": r.trend, "n_years": len(yr),
                 "winter_years": f"{int(yr.index.min())}-{int(yr.index.max())}"})
    print(f"  {lab:6s}: {r.slope*10:+.3f} {unit} ({star(r.p)}, p={r.p:.3f}, n={len(yr)})")
pd.DataFrame(rows).to_csv("TableS1_winter_DecFeb_sensitivity.csv", index=False)

print("\nSaved Table2_annual_trends.csv, TableS1_seasonal_trends.csv and "
      "TableS1_winter_DecFeb_sensitivity.csv")
