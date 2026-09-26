"""
============================================================================
02_merge_and_deseasonalize.py
----------------------------------------------------------------------------
Purpose : Merge the three variables into one master monthly dataset, convert
          the MODIS chlorophyll log-mean to a geometric mean, and compute
          deseasonalized (monthly-anomaly) series for the coupling analysis.
Input   : BoB_EEZ_SST_Chl_monthly_2005_2025.csv   (exported from Google Earth
                                                   Engine; SST + Chl-a)
          BoB_EEZ_SSS_monthly_2005_2025.csv        (from 01_process_SSS.py)
Output  : BoB_EEZ_master_with_anomalies.csv
Notes   : The GEE CSV stores a sentinel value of -999 for fully cloud-covered
          months; these are converted to NaN here.
============================================================================
"""
import numpy as np
import pandas as pd

SC_IN   = "BoB_EEZ_SST_Chl_monthly_2005_2025.csv"   # SST + Chl from GEE
SSS_IN  = "BoB_EEZ_SSS_monthly_2005_2025.csv"       # SSS from script 01
CSV_OUT = "BoB_EEZ_master_with_anomalies.csv"

# ---- 1. SST + Chl (replace -999 sentinel with NaN) ----
sc = pd.read_csv(SC_IN)
for c in ["sst_C", "chl_logmean", "chl_median"]:
    sc[c] = sc[c].replace(-999, np.nan)
# geometric mean chlorophyll = 10 ^ (mean of log10)
sc["chl_geomean"] = 10 ** sc["chl_logmean"]

# ---- 2. SSS ----
ss = pd.read_csv(SSS_IN)

# ---- 3. merge on the year-month key ----
m = sc.merge(ss[["date", "sss_PSU"]], on="date", how="outer").sort_values("date").reset_index(drop=True)

# ---- 4. deseasonalize: anomaly = value - climatological monthly mean ----
for v in ["sst_C", "sss_PSU", "chl_geomean"]:
    clim = m.groupby("month")[v].transform("mean")
    m[v + "_anom"] = m[v] - clim
# chlorophyll anomaly in log10 space (its distribution is log-normal)
m["chl_log"]      = np.log10(m["chl_geomean"])
m["chl_log_anom"] = m["chl_log"] - m.groupby("month")["chl_log"].transform("mean")

cols = ["date", "year", "month", "season",
        "sst_C", "sss_PSU", "chl_geomean", "chl_median", "chl_validpix", "chl_validpct",
        "sst_C_anom", "sss_PSU_anom", "chl_geomean_anom", "chl_log", "chl_log_anom"]
cols = [c for c in cols if c in m.columns]
m[cols].to_csv(CSV_OUT, index=False)

print(f"Saved {CSV_OUT}  ({len(m)} rows)")
print("Missing months per variable:",
      m[["sst_C", "sss_PSU", "chl_geomean"]].isna().sum().to_dict())
print("\nSeasonal climatology:")
print(m.groupby("season")[["sst_C", "sss_PSU", "chl_geomean"]].mean().round(3)
        .reindex(["pre_monsoon", "sw_monsoon", "post_monsoon", "winter"]).to_string())
