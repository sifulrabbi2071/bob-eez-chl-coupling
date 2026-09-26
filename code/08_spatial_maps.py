"""
============================================================================
08_spatial_maps.py
----------------------------------------------------------------------------
Purpose : Seasonal climatology maps and per-pixel Sen's-slope trend maps for
          SST, SSS (GLORYS12V1) and surface Chl-a (Copernicus-GlobColour).
          Field significance of the per-pixel trends is assessed by
          controlling the false discovery rate (Benjamini-Hochberg, q = 0.10;
          Wilks 2016).
Input   : GLORYS_SST_BoB_2005_2025.nc  (thetao)
          GLORYS_SSS_BoB_2005_2025.nc  (so)
          GlobColour_CHL_BoB_2005_2025.nc   (CHL)
          EEZ.shp (+ siblings)
Output  : Fig08-Fig10 climatologies, Fig11a-c trend maps, field_significance_Fig11.csv
Note    : the first run of cartopy downloads coastline/land shapefiles.
============================================================================
"""
import numpy as np
import geopandas as gpd
import xarray as xr
import regionmask
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import LogNorm
from scipy.stats import kendalltau, false_discovery_control
import pandas as pd

mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "savefig.dpi": 300})
EEZ_SHP = "EEZ.shp"
EXT = [88.5, 93.0, 17.5, 23.0]
SEAS = {"pre_monsoon": [3, 4, 5], "sw_monsoon": [6, 7, 8, 9],
        "post_monsoon": [10, 11], "winter": [12, 1, 2]}
SN = {"pre_monsoon": "Pre-monsoon", "sw_monsoon": "SW monsoon",
      "post_monsoon": "Post-monsoon", "winter": "Winter"}
eez = gpd.read_file(EEZ_SHP).to_crs("EPSG:4326")

def load(nc, var):
    ds = xr.open_dataset(nc); da = ds[var]
    for d in ["depth", "deptht", "lev"]:
        if d in da.dims: da = da.isel({d: 0}); break
    da = da.rename({"latitude": "lat", "longitude": "lon"})
    mask = regionmask.mask_geopandas(eez, da["lon"], da["lat"])
    return da.where(~np.isnan(mask))

def _base(ax):
    eez.boundary.plot(ax=ax, edgecolor="k", lw=0.7, transform=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor="#e5e0d8", zorder=2)
    ax.add_feature(cfeature.COASTLINE, lw=0.4)
    ax.set_extent(EXT, crs=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True, lw=0.3, ls=":", alpha=0.5)
    gl.top_labels = False; gl.right_labels = False
    gl.xlabel_style = {"size": 7}; gl.ylabel_style = {"size": 7}

def climatology(da, cmap, label, fname, title, lognorm=False, vmin=None, vmax=None):
    clim = {k: da.sel(time=da["time.month"].isin(mos)).mean("time") for k, mos in SEAS.items()}
    allv = np.concatenate([c.values[~np.isnan(c.values)] for c in clim.values()])
    vmin = np.nanpercentile(allv, 2) if vmin is None else vmin
    vmax = np.nanpercentile(allv, 98) if vmax is None else vmax
    norm = LogNorm(vmin=max(vmin, 0.05), vmax=vmax) if lognorm else None
    fig, axs = plt.subplots(1, 4, figsize=(13, 3.7), subplot_kw={"projection": ccrs.PlateCarree()})
    for ax, (k, mos) in zip(axs, SEAS.items()):
        c = clim[k]
        kw = dict(cmap=cmap, transform=ccrs.PlateCarree(), shading="auto")
        im = (ax.pcolormesh(c["lon"], c["lat"], c.values, norm=norm, **kw) if lognorm
              else ax.pcolormesh(c["lon"], c["lat"], c.values, vmin=vmin, vmax=vmax, **kw))
        _base(ax); ax.set_title(SN[k], fontsize=10, fontweight="bold")
    cb = fig.colorbar(im, ax=axs, orientation="vertical", shrink=0.75, pad=0.02); cb.set_label(label)
    fig.suptitle(title, fontsize=11, fontweight="bold", y=1.02)
    plt.savefig(fname, bbox_inches="tight", facecolor="white"); plt.close(); print("saved", fname)

def sen(y):
    y = np.asarray(y, float); x = np.arange(len(y)); ok = ~np.isnan(y)
    if ok.sum() < 10: return np.nan, np.nan
    xs, ys = x[ok], y[ok]
    sl = [(ys[j]-ys[i])/(xs[j]-xs[i]) for i in range(len(xs)) for j in range(i+1, len(xs))]
    tau, p = kendalltau(xs, ys); return np.median(sl), p

FIELD = []
def trend(da, cmap, label, fname, title, vlim, letter):
    ann = da.groupby("time.year").mean("time")
    lat, lon, arr = ann["lat"].values, ann["lon"].values, ann.values
    slope = np.full((len(lat), len(lon)), np.nan); pval = np.full_like(slope, np.nan)
    for i in range(len(lat)):
        for j in range(len(lon)):
            s, p = sen(arr[:, i, j])
            if not np.isnan(s): slope[i, j] = s*10; pval[i, j] = p
    fig = plt.figure(figsize=(5.4, 5.2)); ax = plt.axes(projection=ccrs.PlateCarree())
    im = ax.pcolormesh(lon, lat, slope, cmap=cmap, vmin=-vlim, vmax=vlim,
                       transform=ccrs.PlateCarree(), shading="auto")
    LO, LA = np.meshgrid(lon, lat); sig = pval < 0.05
    ax.scatter(LO[sig], LA[sig], s=0.6, c="k", alpha=0.35, marker=".", transform=ccrs.PlateCarree())
    _base(ax)
    cb = plt.colorbar(im, ax=ax, shrink=0.8, pad=0.03); cb.set_label(label)
    ax.set_title(title, fontsize=10.5, fontweight="bold")
    fig.text(0.01, 0.99, f"({letter})", fontsize=14, fontweight="bold", va="top")
    ok = np.isfinite(pval); q = np.full_like(pval, np.nan); q[ok] = false_discovery_control(pval[ok])
    FIELD.append({"map": fname, "pixels": int(ok.sum()), "pct_p_lt_0.05": round(100 * (pval[ok] < 0.05).mean(), 1),
                  "pct_FDR_q_lt_0.10": round(100 * (q[ok] < 0.10).mean(), 1)})
    print(f"   field significance: {FIELD[-1]}")
    plt.savefig(fname, bbox_inches="tight", facecolor="white"); plt.close(); print("saved", fname)

# ---- load fields ----
sst = load("GLORYS_SST_BoB_2005_2025.nc", "thetao")
sss = load("GLORYS_SSS_BoB_2005_2025.nc", "so")
chl = load("GlobColour_CHL_BoB_2005_2025.nc", "CHL")

# ---- climatologies ----
climatology(sst, "inferno", "SST (°C)", "Fig08_SST_climatology.png",
            "Seasonal climatology of sea surface temperature (2005-2025)")
climatology(sss, "viridis_r", "SSS (PSU)", "Fig09_SSS_climatology.png",
            "Seasonal climatology of sea surface salinity (2005-2025)", vmin=19, vmax=32)
climatology(chl, "YlGn", "Chl-a (mg m$^{-3}$)", "Fig10_CHL_climatology.png",
            "Seasonal climatology of surface chlorophyll-a (2005-2025)", lognorm=True, vmin=0.1, vmax=5)

# ---- trends ----
trend(sst, "RdBu_r", "SST trend (°C decade$^{-1}$)", "Fig11a_SST_trend.png",
      "Sea surface temperature trend (2005-2025)", vlim=0.6, letter="a")
trend(sss, "RdBu", "SSS trend (PSU decade$^{-1}$)", "Fig11b_SSS_trend.png",
      "Sea surface salinity trend (2005-2025)", vlim=1.0, letter="b")
trend(np.log10(chl.where(chl > 0)), "BrBG", "log$_{10}$(Chl-a) trend (decade$^{-1}$)",
      "Fig11c_CHL_trend.png", "Surface chlorophyll-a trend (2005-2025)", vlim=0.15, letter="c")

pd.DataFrame(FIELD).to_csv("field_significance_Fig11.csv", index=False)
print("Saved field_significance_Fig11.csv")
print("All maps done.")
