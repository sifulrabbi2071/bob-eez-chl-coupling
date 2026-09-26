"""
============================================================================
01_process_SSS.py
----------------------------------------------------------------------------
Purpose : Convert the GLORYS12V1 sea-surface-salinity NetCDF into a monthly,
          EEZ-averaged time series (2005-2025).
Input   : GLORYS_SSS_BoB_2005_2025.nc   (downloaded from Copernicus Marine)
          EEZ.shp (+ .shx .dbf .prj)     (Bangladesh EEZ polygon; created by 00b_get_EEZ.py)
Output  : BoB_EEZ_SSS_monthly_2005_2025.csv
============================================================================
"""
import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import regionmask

# ---- paths (edit if your files live elsewhere) ----
NC_IN    = "GLORYS_SSS_BoB_2005_2025.nc"
EEZ_SHP  = "EEZ.shp"
CSV_OUT  = "BoB_EEZ_SSS_monthly_2005_2025.csv"

# ---- 1. load salinity, keep the surface level ----
ds  = xr.open_dataset(NC_IN)
sss = ds["so"]                        # sea_water_salinity, PSU
for d in ["depth", "deptht", "lev"]:  # drop the (single) depth level if present
    if d in sss.dims:
        sss = sss.isel({d: 0})
        break
sss = sss.rename({"latitude": "lat", "longitude": "lon"})

# ---- 2. mask to the Bangladesh EEZ polygon ----
eez  = gpd.read_file(EEZ_SHP).to_crs("EPSG:4326")
mask = regionmask.mask_geopandas(eez, sss["lon"], sss["lat"])   # NaN outside EEZ
sss_eez = sss.where(~np.isnan(mask))

# ---- 3. cosine-latitude-weighted monthly EEZ mean ----
weights   = np.cos(np.deg2rad(sss_eez["lat"]))
sss_series = sss_eez.weighted(weights).mean(dim=["lat", "lon"])

# ---- 4. tidy into a dataframe with a monsoon-season label ----
df = sss_series.to_dataframe(name="sss_PSU").reset_index()
df["time"]  = pd.to_datetime(df["time"])
df["year"]  = df["time"].dt.year
df["month"] = df["time"].dt.month
df["date"]  = df["time"].dt.strftime("%Y-%m")

def bob_season(m):
    if m in (3, 4, 5):    return "pre_monsoon"
    if m in (6, 7, 8, 9): return "sw_monsoon"
    if m in (10, 11):     return "post_monsoon"
    return "winter"                       # Dec, Jan, Feb
df["season"] = df["month"].apply(bob_season)

df = df[["date", "year", "month", "season", "sss_PSU"]].sort_values("date").reset_index(drop=True)
df.to_csv(CSV_OUT, index=False)

print(f"Saved {CSV_OUT}  ({len(df)} rows)")
print("SSS (PSU) summary:", df["sss_PSU"].describe()[["min", "mean", "max"]].round(2).to_dict())
print("\nSeasonal mean SSS:")
print(df.groupby("season")["sss_PSU"].mean().round(2).to_string())
