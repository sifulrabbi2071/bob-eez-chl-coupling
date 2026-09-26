"""
============================================================================
05_climate_indices.py
----------------------------------------------------------------------------
Purpose : Relate the deseasonalized anomalies to the Indian Ocean Dipole
          (DMI) and ENSO (Nino 3.4): contemporaneous, lagged, and seasonal
          correlations. Because the indices are strongly autocorrelated,
          lagged correlations are tested with an effective sample size
          (Bretherton et al. 1999) and seasonal correlations with the 21
          seasonal means (one value per year).
Input   : BoB_EEZ_master_with_anomalies.csv   (from script 02)
          Climate indices are downloaded automatically from NOAA PSL.
Output  : climate_indices_2005_2025.csv, BoB_EEZ_master_with_indices.csv,
          climate_mode_correlations.csv  + console summary
============================================================================
"""
import urllib.request
import numpy as np
import pandas as pd
from scipy import stats

M_IN  = "BoB_EEZ_master_with_anomalies.csv"
DMI_URL   = "https://psl.noaa.gov/gcos_wgsp/Timeseries/Data/dmi.had.long.data"
NINO_URL  = "https://psl.noaa.gov/data/correlation/nina34.anom.data"

def fetch(url):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    return urllib.request.urlopen(req, timeout=30).read().decode("utf-8", "ignore")

def parse_psl(txt):
    """Parse a NOAA PSL year + 12-monthly-value table."""
    lines = txt.splitlines()
    y0, y1 = (int(x) for x in lines[0].split()[:2])
    rows = []
    for ln in lines[1:]:
        p = ln.split()
        if len(p) >= 13:
            try:
                yr = int(p[0]); vals = [float(x) for x in p[1:13]]
            except ValueError:
                continue
            if yr < y0 or yr > y1:
                continue
            for mth, v in enumerate(vals, 1):
                rows.append((yr, mth, v))
    d = pd.DataFrame(rows, columns=["year", "month", "val"])
    d.loc[(d["val"] <= -9) | (d["val"] >= 99), "val"] = np.nan   # missing flags
    return d

def eff_p(a, b):
    """Pearson r and p using an effective sample size from lag-1 autocorrelations."""
    q = pd.DataFrame({"a": a, "b": b}).dropna()
    r = q["a"].corr(q["b"]); r1, r2 = q["a"].autocorr(1), q["b"].autocorr(1)
    n_eff = len(q) * (1 - r1 * r2) / (1 + r1 * r2)
    t = r * np.sqrt((n_eff - 2) / (1 - r**2))
    return r, 2 * stats.t.sf(abs(t), n_eff - 2), n_eff

def star(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"

# ---- 1. download + assemble indices ----
dmi  = parse_psl(fetch(DMI_URL)).rename(columns={"val": "IOD"})
nino = parse_psl(fetch(NINO_URL)).rename(columns={"val": "ENSO"})
idx = dmi.merge(nino, on=["year", "month"], how="outer")
idx = idx[(idx["year"] >= 2005) & (idx["year"] <= 2025)].sort_values(["year", "month"]).reset_index(drop=True)
idx.to_csv("climate_indices_2005_2025.csv", index=False)

# ---- 2. merge with the ocean anomalies ----
m = pd.read_csv(M_IN).merge(idx[["year", "month", "IOD", "ENSO"]], on=["year", "month"], how="left")
m = m.sort_values("date").reset_index(drop=True)
m.to_csv("BoB_EEZ_master_with_indices.csv", index=False)

avars = {"SST": "sst_C_anom", "SSS": "sss_PSU_anom", "Chl": "chl_log_anom"}

print("=" * 70)
print("A. CONTEMPORANEOUS correlation (anomaly vs index)")
print("=" * 70)
for name, col in avars.items():
    d = m.dropna(subset=[col, "IOD", "ENSO"])
    r1, p1 = stats.pearsonr(d[col], d["IOD"])
    r2, p2 = stats.pearsonr(d[col], d["ENSO"])
    print(f"  {name:4s}: IOD {r1:+.2f}({star(p1)})   ENSO {r2:+.2f}({star(p2)})")

print("\n" + "=" * 70)
print("B. LAGGED correlation (index leads by k months; p from effective sample size)")
print("=" * 70)
res = []
print(f"{'lag':4s}{'IOD-Chl':>12s}{'ENSO-Chl':>12s}{'IOD-SSS':>12s}")
for lag in range(0, 8):
    out = []
    for ic, vc in [("IOD", "chl_log_anom"), ("ENSO", "chl_log_anom"), ("IOD", "sss_PSU_anom")]:
        r, p, ne = eff_p(m[ic].shift(lag).values, m[vc].values)
        out.append(f"{r:+.2f}{star(p)}")
        res.append({"analysis": "lagged", "lag_or_season": lag, "pair": f"{ic}-{vc}", "r": round(r, 3),
                    "p": round(p, 4), "n_eff": round(ne, 1)})
    print(f"{lag:<4d}{out[0]:>12s}{out[1]:>12s}{out[2]:>12s}")

print("\n" + "=" * 70)
print("C. SEASONAL correlation (r from monthly anomalies; p from the 21 seasonal means)")
print("=" * 70)
for s in ["pre_monsoon", "sw_monsoon", "post_monsoon", "winter"]:
    d = m[m["season"] == s]
    out = []
    for ic, vc in [("IOD", "chl_log_anom"), ("ENSO", "chl_log_anom"),
                   ("IOD", "sss_PSU_anom"), ("ENSO", "sss_PSU_anom")]:
        dd = d.dropna(subset=[ic, vc])
        r = stats.pearsonr(dd[ic], dd[vc])[0]
        y = dd.groupby("year")[[ic, vc]].mean().dropna()
        ry, p = stats.pearsonr(y[ic], y[vc])
        out.append(f"{r:+.2f}{star(p)}")
        res.append({"analysis": "seasonal", "lag_or_season": s, "pair": f"{ic}-{vc}", "r": round(r, 3),
                    "r_seasonal_means": round(ry, 3), "p": round(p, 4), "n_eff": len(y)})
    print(f"  {s:13s}: IOD-Chl {out[0]}  ENSO-Chl {out[1]}  IOD-SSS {out[2]}  ENSO-SSS {out[3]}")

pd.DataFrame(res).to_csv("climate_mode_correlations.csv", index=False)
print("\nSaved climate_indices_2005_2025.csv, BoB_EEZ_master_with_indices.csv, climate_mode_correlations.csv")
