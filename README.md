# Seasonally shifting thermal–haline control of surface chlorophyll-a in the Bangladesh EEZ

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Analysis code and derived data for the study of 21 years (2005–2025) of sea
surface temperature (SST), sea surface salinity (SSS) and surface chlorophyll-a
(Chl-a) in the Bangladesh Exclusive Economic Zone (EEZ), northern Bay of Bengal.

> **Before publishing:** add the GitHub URL, release date and Zenodo DOI to
> `CITATION.cff`, and the paper DOI below once available.

---

## Overview

The workflow combines **Google Earth Engine** (JavaScript; monthly EEZ-mean
SST and MODIS Chl-a) with **Python** (all other processing, statistics and figures). Running
the scripts in order reproduces every number, table and figure in the paper.

**Main results reproduced by this code**
- EEZ warming of +0.32 °C decade⁻¹ (p < 0.001) and freshening of −0.36 PSU decade⁻¹
  (p = 0.027); no significant Chl-a trend (+0.04 mg m⁻³ decade⁻¹, p = 0.10).
  All from the Hirsch–Slack seasonal Mann–Kendall test.
- A seasonal handover in the dominant partial association with surface Chl-a:
  salinity in the SW monsoon and post-monsoon, temperature in winter.
- Robustness of the handover to cloud-related data gaps (coverage thresholds of
  10%, 15% and 20%), to serial dependence (seasonal-mean test) and to excluding
  turbid coastal water (offshore test).
- Fresher post-monsoon water during positive Indian Ocean Dipole phases; weak
  links between the climate modes and Chl-a.
- A spatially coherent dipole in the Chl-a trend that remains significant after
  false-discovery-rate correction: increases along the coast, decreases offshore.

---

## Repository structure

```
.
├── code/                        analysis scripts (run in the order below)
│   ├── 00b_get_EEZ.py                   EEZ + Bay of Bengal polygons (Marine Regions)
│   ├── 00_GEE_SST_Chl_extraction.js     Google Earth Engine: SST + MODIS Chl-a series
│   ├── 01_process_SSS.py                GLORYS SSS -> EEZ-mean series
│   ├── 02_merge_and_deseasonalize.py    master dataset + anomalies
│   ├── 03_trend_analysis.py             Table 2, Table S1, Dec–Feb sensitivity
│   ├── 03b_seasonal_autocorrelation.py  Table S1b (autocorrelation of seasonal series)
│   ├── 04_coupling_analysis.py          Table 3 / S2 (10/15/20% thresholds), S3, S4a
│   ├── 04b_offshore_test.py             Table S4b (offshore sub-domain test)
│   ├── 05_climate_indices.py            IOD / ENSO correlations
│   ├── 06_globcolour_crosscheck.py      MODIS vs Copernicus-GlobColour
│   ├── 07_figures_main.py               Figs. 3–7
│   ├── 08_spatial_maps.py               Figs. 8–11 + field significance
│   ├── 09_study_area_and_workflow.py    Figs. 1–2
│   └── 10_supporting_numbers.py         other numbers quoted in the text
├── data/        derived monthly EEZ time series
├── tables/      Table 2 and Supplementary Tables S1, S1b, S2–S4 (CSV)
├── figures/     all figures, numbered as in the paper
├── requirements.txt, LICENSE, CITATION.cff, CHANGELOG.md
└── README.md
```

---

## Data

Raw source data are **not redistributed**; each product is freely available from
its provider. The `data/` folder holds the **derived** EEZ-mean monthly series,
so the statistics and EEZ-mean figures can be regenerated without the raw files.

| Variable | Product | Provider |
|----------|---------|----------|
| SST (EEZ series) | NOAA OISST v2.1 | NOAA NCEI, via Google Earth Engine |
| Chl-a (EEZ series) | MODIS-Aqua L3SMI (R2022.0) | NASA OBPG, via Google Earth Engine |
| SSS; SST for maps | GLORYS12V1 reanalysis | Copernicus Marine Service |
| Chl-a (cross-check, maps) | Copernicus-GlobColour L4 (OCEANCOLOUR_GLO_BGC_L4_MY_009_104) | Copernicus Marine Service |
| EEZ polygon | Maritime Boundaries Geodatabase v12, MRGID 8481 | Flanders Marine Institute (marineregions.org) |
| IOD, ENSO indices | DMI (HadISST), Niño 3.4 | NOAA PSL (downloaded by script 05) |

### Copernicus Marine downloads (free account required)
```bash
copernicusmarine subset --dataset-id cmems_mod_glo_phy_my_0.083deg_P1M-m --variable so     --minimum-longitude 87.5 --maximum-longitude 93.0 --minimum-latitude 15.0 --maximum-latitude 23.0 --minimum-depth 0 --maximum-depth 1 --start-datetime 2005-01-01 --end-datetime 2025-12-31 --output-filename GLORYS_SSS_BoB_2005_2025.nc --output-directory .
copernicusmarine subset --dataset-id cmems_mod_glo_phy_my_0.083deg_P1M-m --variable thetao --minimum-longitude 87.5 --maximum-longitude 93.0 --minimum-latitude 15.0 --maximum-latitude 23.0 --minimum-depth 0 --maximum-depth 1 --start-datetime 2005-01-01 --end-datetime 2025-12-31 --output-filename GLORYS_SST_BoB_2005_2025.nc --output-directory .
copernicusmarine subset --dataset-id cmems_obs-oc_glo_bgc-plankton_my_l4-multi-4km_P1M --variable CHL --minimum-longitude 87.5 --maximum-longitude 93.0 --minimum-latitude 15.0 --maximum-latitude 23.0 --start-datetime 2005-01-01 --end-datetime 2025-12-31 --output-filename GlobColour_CHL_BoB_2005_2025.nc --output-directory .
```
Log in with `copernicusmarine login`. Never store your credentials in the repository.

---

## How to reproduce

All scripts read and write files in the **current working folder**. Run them
from a folder containing the downloaded NetCDF files.

```bash
# 1. environment
conda create -n boban -c conda-forge python=3.12 --file requirements.txt
conda activate boban

# 2. boundary polygons -> EEZ.shp, EEZ.zip, bob_iho.geojson
python code/00b_get_EEZ.py

# 3. Google Earth Engine (browser):
#    upload EEZ.zip as a table asset in your own Earth Engine project, put its
#    ID in EEZ_ASSET at the top of code/00_GEE_SST_Chl_extraction.js (the
#    authors' asset is used by default), run the script, and copy the exported
#    BoB_EEZ_SST_Chl_monthly_2005_2025.csv into the working folder
#    (or use the copy in data/).

# 4. Python pipeline
python code/01_process_SSS.py
python code/02_merge_and_deseasonalize.py
python code/03_trend_analysis.py
python code/03b_seasonal_autocorrelation.py
python code/04_coupling_analysis.py
python code/04b_offshore_test.py
python code/05_climate_indices.py
python code/06_globcolour_crosscheck.py
python code/07_figures_main.py
python code/08_spatial_maps.py
python code/09_study_area_and_workflow.py
python code/10_supporting_numbers.py
```
Scripts 02–05 (including 03b) and 07 need only the CSV files in `data/`. Scripts 01, 04b, 06,
08, 09 and 10 also need the NetCDF files and/or `EEZ.shp`.

---

## Where each result comes from

| Paper item | Script | Output |
|------------|--------|--------|
| Table 1 | — | product list (this README) |
| Table 2 | 03 | `tables/Table2_annual_trends.csv` |
| Table 3, Fig. 5 | 04, 07 | `tables/TableS2_seasonal_coupling.csv` (full-record columns) |
| Fig. 1, Fig. 2 | 09 | `figures/Fig01_*`, `Fig02_*` |
| Fig. 3, Fig. 4 | 03, 07 | `figures/Fig03_*`, `Fig04_*`; `tables/TableS1_*` |
| Fig. 6 | 05, 07 | `figures/Fig06_*`; `tables/climate_mode_correlations.csv` |
| Fig. 7 (§3.6) | 06, 07 | `figures/Fig07_*` |
| Figs. 8–11 (§3.7) | 08 | `figures/Fig08_*`–`Fig11c_*`; `tables/field_significance_Fig11.csv` |
| Table S1 | 03 | seasonal trends; continuous Dec–Feb winter sensitivity (`tables/TableS1_winter_DecFeb_sensitivity.csv`, §3.2) |
| Table S1b (§2.5) | 03b | lag-1 autocorrelation of the detrended seasonal series |
| Tables S2, S3 (§2.6, §3.4) | 04 | coverage-threshold test (10%, 15%, 20%); Chl-a coverage by season |
| Table S4 (a, b) | 04, 04b | seasonal-mean test; offshore test (§3.4) |
| Other numbers in §2.1, §2.3, §2.7, §3.1, §3.7 | 10 | `tables/supporting_numbers.txt` |

---

## Software

Python 3.12 with NumPy, pandas, SciPy, xarray, netCDF4, geopandas, regionmask,
pyMannKendall, Matplotlib and Cartopy. Tested versions are listed in
`requirements.txt`.

## Citation

Please cite the paper and this repository (see `CITATION.cff`).

## License

Code and derived data are released under the MIT License (see `LICENSE`). The
source datasets remain subject to their providers' licensing terms.

## Author and contact

**Md Siful Islam Rabbi**, Noakhali Science and Technology University, Bangladesh
Email: sifulrabbi2871@gmail.com
