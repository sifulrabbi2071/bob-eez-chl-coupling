"""
============================================================================
04b_offshore_test.py
----------------------------------------------------------------------------
Purpose : Test whether the seasonal couplings could arise from optically
          complex coastal water (CDOM, suspended sediment). The partial
          correlations are recomputed with a product combination available on
          a grid (GLORYS12V1 SST and SSS, Copernicus-GlobColour Chl-a), first
          for the full EEZ and then for the EEZ south of 21 N.
Input   : GLORYS_SST_BoB_2005_2025.nc, GLORYS_SSS_BoB_2005_2025.nc,
          GlobColour_CHL_BoB_2005_2025.nc, EEZ.shp,
          BoB_EEZ_master_with_anomalies.csv  (for the season labels)
Output  : TableS4b_coupling_offshore.csv  + console summary
============================================================================
"""
import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import regionmask
from scipy import stats
from shapely.geometry import box

LAT_CUT = 21.0
eez = gpd.read_file("EEZ.shp").to_crs("EPSG:4326")
lab = pd.read_csv("BoB_EEZ_master_with_anomalies.csv")[["date", "year", "month", "season"]]
SEASONS = ["pre_monsoon", "sw_monsoon", "post_monsoon", "winter"]

def star(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"

def load(nc, var):
    da = xr.open_dataset(nc)[var]
    for d in ["depth", "deptht", "lev"]:
        if d in da.dims: da = da.isel({d: 0}); break
    da = da.rename({"latitude": "lat", "longitude": "lon"})
    return da.where(~np.isnan(regionmask.mask_geopandas(eez, da["lon"], da["lat"])))

def eez_mean(da, lat_max, log=False):
    da = da.where(da["lat"] < lat_max)
    if log: da = np.log10(da.where(da > 0))
    return da.weighted(np.cos(np.deg2rad(da["lat"]))).mean(["lat", "lon"]).values

def partial_corr(x, y, z):
    def resid(a, b):
        b = np.column_stack([np.ones_like(b), b]); return a - b @ np.linalg.lstsq(b, a, rcond=None)[0]
    r = np.corrcoef(resid(x, z), resid(y, z))[0, 1]; dof = len(x) - 3
    return r, 2 * stats.t.sf(abs(r * np.sqrt(dof / (1 - r**2))), dof)

sst = load("GLORYS_SST_BoB_2005_2025.nc", "thetao")
sss = load("GLORYS_SSS_BoB_2005_2025.nc", "so")
chl = load("GlobColour_CHL_BoB_2005_2025.nc", "CHL")
aea = "+proj=aea +lat_1=15 +lat_2=23 +lon_0=90"                 # equal-area projection
frac = float(eez.clip(box(80, 0, 100, LAT_CUT)).to_crs(aea).area.sum() / eez.to_crs(aea).area.sum())
print(f"Offshore sub-domain (south of {LAT_CUT} N) = {frac*100:.0f}% of the EEZ area")

rows = []
for dom, lat_max in [("full EEZ", 99.0), (f"south of {LAT_CUT:g} N", LAT_CUT)]:
    d = lab.copy()
    d["T"], d["S"], d["C"] = eez_mean(sst, lat_max), eez_mean(sss, lat_max), eez_mean(chl, lat_max, log=True)
    for v in "TSC":
        d[v + "a"] = d[v] - d.groupby("month")[v].transform("mean")
    print(f"\n{dom}")
    for s in SEASONS:
        q = d[d["season"] == s].dropna(subset=["Ta", "Sa", "Ca"])
        for a, ctrl, name in [("Ta", "Sa", "SST-Chl (control SSS)"), ("Sa", "Ta", "SSS-Chl (control SST)")]:
            r, p = partial_corr(q[a].values, q["Ca"].values, q[ctrl].values)
            rows.append({"Domain": dom, "Season": s, "Coupling": name, "partial_r": round(r, 2),
                         "p": round(p, 4), "sig": star(p), "n": len(q)})
            print(f"  {s:13s} {name:22s}: {r:+.2f} ({star(p)}, p={p:.3f})")
pd.DataFrame(rows).to_csv("TableS4b_coupling_offshore.csv", index=False)
print("\nSaved TableS4b_coupling_offshore.csv")
