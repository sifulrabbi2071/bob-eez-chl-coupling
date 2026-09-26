"""
============================================================================
06_globcolour_crosscheck.py
----------------------------------------------------------------------------
Purpose : Cross-validate the MODIS-Aqua chlorophyll trend against the merged,
          gap-filled Copernicus-GlobColour L4 product (CMEMS
          OCEANCOLOUR_GLO_BGC_L4_MY_009_104), processed with the same pipeline.
          Note: GlobColour includes MODIS-Aqua among its input sensors, so the
          comparison is only partly independent.
Input   : GlobColour_CHL_BoB_2005_2025.nc      (CMEMS, monthly, 4 km)
          EEZ.shp (+ siblings)
          BoB_EEZ_master_with_anomalies.csv    (for the MODIS series)
Output  : BoB_EEZ_GlobColour_monthly_2005_2025.csv
          BoB_EEZ_chl_MODIS_vs_GlobColour.csv  + console summary
============================================================================
"""
import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import regionmask
import pymannkendall as mk

NC_IN   = "GlobColour_CHL_BoB_2005_2025.nc"
EEZ_SHP = "EEZ.shp"
MOD_IN  = "BoB_EEZ_master_with_anomalies.csv"

def star(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"

# ---- 1. GlobColour: EEZ clip + geometric-mean averaging (matches MODIS pipeline) ----
chl = xr.open_dataset(NC_IN)["CHL"].rename({"latitude": "lat", "longitude": "lon"})
eez = gpd.read_file(EEZ_SHP).to_crs("EPSG:4326")
chl_eez = chl.where(~np.isnan(regionmask.mask_geopandas(eez, chl["lon"], chl["lat"])))
log_chl = np.log10(chl_eez.where(chl_eez > 0))
geomean = 10 ** log_chl.weighted(np.cos(np.deg2rad(chl_eez["lat"]))).mean(dim=["lat", "lon"])

df = geomean.to_dataframe(name="chl_globcolour").reset_index()
df["time"] = pd.to_datetime(df["time"])
df["year"], df["month"] = df["time"].dt.year, df["time"].dt.month
df["date"] = df["time"].dt.strftime("%Y-%m")
df["season"] = df["month"].map(lambda k: "pre_monsoon" if k in (3, 4, 5) else "sw_monsoon" if k in (6, 7, 8, 9)
                               else "post_monsoon" if k in (10, 11) else "winter")
df = df[["date", "year", "month", "season", "chl_globcolour"]].sort_values("date").reset_index(drop=True)
df.to_csv("BoB_EEZ_GlobColour_monthly_2005_2025.csv", index=False)

# ---- 2. MODIS series ----
mod = pd.read_csv(MOD_IN)[["date", "year", "month", "season", "chl_geomean"]].rename(columns={"chl_geomean": "chl_modis"})
both = mod.merge(df[["date", "chl_globcolour"]], on="date", how="outer").sort_values("date").reset_index(drop=True)
both.to_csv("BoB_EEZ_chl_MODIS_vs_GlobColour.csv", index=False)

# ---- 3. compare full-record trends (Hirsch-Slack; NaN months kept for calendar alignment) ----
print("=" * 70)
print("FULL-RECORD CHLOROPHYLL TREND  (Hirsch-Slack seasonal Mann-Kendall)")
print("=" * 70)
for name, col in [("MODIS-Aqua", "chl_modis"), ("GlobColour", "chl_globcolour")]:
    r = mk.correlated_seasonal_test(both[col].values, period=12)
    print(f"  {name:11s}: {r.slope*10:+.4f} mg/m3/decade   p={r.p:.4f}   {r.trend}")

c = both.dropna(subset=["chl_modis", "chl_globcolour"])
print(f"\nMonthly correlation MODIS vs GlobColour: r = {np.corrcoef(c['chl_modis'], c['chl_globcolour'])[0,1]:.3f} (n={len(c)})")
print(f"GlobColour mean / MODIS mean = {c['chl_globcolour'].mean()/c['chl_modis'].mean():.2f}")
print("\nSeasonal trends (MK on yearly seasonal means):")
for s in ["pre_monsoon", "sw_monsoon", "post_monsoon", "winter"]:
    out = []
    for col in ["chl_modis", "chl_globcolour"]:
        yr = both[both["season"] == s].dropna(subset=[col]).groupby("year")[col].mean()
        r = mk.original_test(yr.values); out.append(f"{r.slope*10:+.3f} ({star(r.p)})")
    print(f"  {s:13s}  MODIS {out[0]:14s} GlobColour {out[1]}")
print("\nSaved BoB_EEZ_GlobColour_monthly_2005_2025.csv and BoB_EEZ_chl_MODIS_vs_GlobColour.csv")
