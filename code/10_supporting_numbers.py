"""
============================================================================
10_supporting_numbers.py
----------------------------------------------------------------------------
Purpose : Supporting numbers quoted in the text that are not part of a table:
            §2.1  share of valid ocean-colour cells in the estuarine channels
            §2.3  effect of cosine-latitude weighting on EEZ means
            §2.7  GLORYS12V1 EEZ-mean SST trend (vs OISST)
            §3.1  record means and ranges
            §3.7  Chl-a trend dipole by zone, before and after FDR control
Input   : BoB_EEZ_master_with_anomalies.csv, EEZ.shp,
          GLORYS_SST_BoB_2005_2025.nc, GlobColour_CHL_BoB_2005_2025.nc
Output  : supporting_numbers.txt  (also printed)
============================================================================
"""
import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import regionmask
import pymannkendall as mk
from scipy.stats import kendalltau, false_discovery_control

out = []
def say(s): print(s); out.append(s)
eez = gpd.read_file("EEZ.shp").to_crs("EPSG:4326")

def load(nc, var):
    da = xr.open_dataset(nc)[var]
    for d in ["depth", "deptht", "lev"]:
        if d in da.dims: da = da.isel({d: 0}); break
    da = da.rename({"latitude": "lat", "longitude": "lon"})
    return da.where(~np.isnan(regionmask.mask_geopandas(eez, da["lon"], da["lat"])))

sst = load("GLORYS_SST_BoB_2005_2025.nc", "thetao")
chl = load("GlobColour_CHL_BoB_2005_2025.nc", "CHL")

# §3.1 record statistics
m = pd.read_csv("BoB_EEZ_master_with_anomalies.csv")
say(f"§3.1  SST mean {m.sst_C.mean():.1f} °C (range {m.sst_C.min():.1f}–{m.sst_C.max():.1f}); "
    f"SSS mean {m.sss_PSU.mean():.1f} (range {m.sss_PSU.min():.1f}–{m.sss_PSU.max():.1f}); "
    f"Chl-a mean {m.chl_geomean.mean():.2f} mg m-3; valid Chl-a months {m.chl_geomean.notna().sum()}")

# §2.1 estuarine channels (north of 22.3 N) as a share of valid ocean-colour cells
valid = chl.isel(time=0).notnull()
say(f"§2.1  estuarine-channel cells: {100 * float((valid & (chl['lat'] > 22.3)).sum() / valid.sum()):.1f}% of valid ocean-colour cells")

# §2.3 weighted vs unweighted EEZ means
def wdiff(da, log=False):
    x = np.log10(da.where(da > 0)) if log else da
    a = x.weighted(np.cos(np.deg2rad(x["lat"]))).mean(["lat", "lon"]); b = x.mean(["lat", "lon"])
    return float(np.nanmax(np.abs((a - b).values)))
say(f"§2.3  max |weighted - unweighted| EEZ mean: SST {wdiff(sst):.4f} °C; log10 Chl-a {wdiff(chl, True):.4f}")

# §2.7 GLORYS EEZ-mean SST trend
g = sst.weighted(np.cos(np.deg2rad(sst["lat"]))).mean(["lat", "lon"]).values
r = mk.correlated_seasonal_test(g, period=12)
say(f"§2.7  GLORYS12V1 EEZ-mean SST trend {r.slope*10:+.2f} °C/decade (Hirsch-Slack p={r.p:.4f})")

# §3.7 Chl-a trend dipole by zone (same method as 08_spatial_maps.py)
ann = np.log10(chl.where(chl > 0)).groupby("time.year").mean("time")
arr = ann.values; LA = np.broadcast_to(ann["lat"].values[:, None], arr.shape[1:])
P, S, L = [], [], []
for i, j in zip(*np.where(np.isfinite(arr).sum(0) >= 10)):
    y = arr[:, i, j]; ok = np.isfinite(y); x = np.arange(len(y))[ok]; y = y[ok]
    S.append(np.median([(y[b] - y[a]) / (x[b] - x[a]) for a in range(len(x)) for b in range(a + 1, len(x))]))
    P.append(kendalltau(x, y)[1]); L.append(LA[i, j])
P, S, L = map(np.array, (P, S, L)); fdr = false_discovery_control(P) < 0.10; raw = P < 0.05
north, south = L >= 21, L < 20
say(f"§3.7  Chl-a significant increases north of 21 N: {100*(raw&(S>0)&north).sum()/north.sum():.0f}% "
    f"(after FDR {100*(fdr&(S>0)&north).sum()/north.sum():.0f}%); significant decreases south of 20 N: "
    f"{100*(raw&(S<0)&south).sum()/south.sum():.0f}% (after FDR {100*(fdr&(S<0)&south).sum()/south.sum():.0f}%)")

open("supporting_numbers.txt", "w").write("\n".join(out) + "\n")
