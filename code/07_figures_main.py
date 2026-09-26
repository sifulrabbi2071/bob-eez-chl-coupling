"""
============================================================================
07_figures_main.py
----------------------------------------------------------------------------
Purpose : EEZ-averaged figures, numbered as in the paper:
            Fig03_timeseries_trend.png    monthly series + seasonal Sen's slope
                                          (Hirsch-Slack statistics, as Table 2)
            Fig04_seasonal_trend.png      seasonal trend bars
            Fig05_coupling_heatmap.png    seasonal partial-correlation map
            Fig06_climate_modes.png       IOD / ENSO lagged + seasonal correlation
            Fig07_crossvalidation.png     MODIS-Aqua vs Copernicus-GlobColour
Input   : BoB_EEZ_master_with_anomalies.csv, BoB_EEZ_master_with_indices.csv,
          BoB_EEZ_chl_MODIS_vs_GlobColour.csv
============================================================================
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy import stats
import pymannkendall as mk

m = pd.read_csv("BoB_EEZ_master_with_anomalies.csv").sort_values("date").reset_index(drop=True)
m["t"] = pd.to_datetime(m["date"])
SEAS  = ["pre_monsoon", "sw_monsoon", "post_monsoon", "winter"]
SNAME = {"pre_monsoon": "Pre-monsoon", "sw_monsoon": "SW monsoon", "post_monsoon": "Post-monsoon", "winter": "Winter"}
def star(p): return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
def minus(s): return s.replace("-0", "\u22120").replace("-1", "\u22121")

# ---------- Fig. 3: time series + Hirsch-Slack seasonal Sen's slope ----------
panels = [("sst_C", "SST (°C)", "#d6604d", "°C", False), ("sss_PSU", "SSS (PSU)", "#1b8a8f", "PSU", False),
          ("chl_geomean", "Chl-a (mg m$^{-3}$)", "#4a9c4e", "mg m$^{-3}$", True)]
fig, axes = plt.subplots(3, 1, figsize=(9, 8.1), sharex=True)
idx = np.arange(len(m))
for ax, (col, ylab, colr, unit, logy) in zip(axes, panels):
    x = m[col].values.astype(float)
    r = mk.correlated_seasonal_test(x, period=12)
    valid = ~np.isnan(x)
    line = (np.nanmedian(x) - np.median(idx[valid]) / 12 * r.slope) + r.slope * idx / 12
    ptxt = "p < 0.001" if r.p < 0.001 else f"p = {r.p:.3f}"
    ax.plot(m["t"], x, color=colr, lw=1)
    ax.plot(m["t"], line, "k--", lw=1.6, label=minus(f"Sen's slope = {r.slope*10:+.2f} {unit} decade$^{{-1}}$ (seasonal MK, {ptxt})"))
    ax.set_ylabel(ylab)
    if logy: ax.set_yscale("log"); ax.set_ylim(0.1, 8)
    else:
        lo, hi = np.nanmin(x), np.nanmax(x); ax.set_ylim(lo - 0.05 * (hi - lo), hi + 0.22 * (hi - lo))
    ax.legend(loc="upper left", frameon=False, fontsize=9); ax.grid(ls=":", alpha=0.5)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
axes[-1].set_xlim(pd.Timestamp("2005-01-01"), pd.Timestamp("2025-12-31"))
axes[-1].xaxis.set_major_locator(mdates.YearLocator(5, month=1, day=1))
axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%Y")); axes[-1].set_xlabel("Year")
fig.suptitle("Monthly time series with Sen's slope trend (2005–2025)", fontweight="bold")
fig.tight_layout(); fig.savefig("Fig03_timeseries_trend.png", dpi=400); plt.close()

# ---------- Figs. 4-5: original plotting style ----------
STYLE = {"font.family": "DejaVu Sans", "font.size": 10, "axes.linewidth": 0.8, "axes.spines.top": False,
         "axes.spines.right": False, "savefig.dpi": 300, "axes.titleweight": "bold"}
C = {"sst": "#d1495b", "sss": "#00798c", "chl": "#2e933c"}
with plt.rc_context(STYLE):
    rows = []
    for v, lab in [("sst_C", "SST"), ("sss_PSU", "SSS"), ("chl_geomean", "Chl-a")]:
        for s in SEAS:
            yr = m[m["season"] == s].dropna(subset=[v]).groupby("year")[v].mean()
            r = mk.original_test(yr.values); rows.append((lab, s, r.slope * 10, r.p))
    td = pd.DataFrame(rows, columns=["var", "season", "slope", "p"])
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.8))
    for ax, (v, col, unit) in zip(axs, [("SST", C["sst"], "°C/dec"), ("SSS", C["sss"], "PSU/dec"),
                                        ("Chl-a", C["chl"], "mg m$^{-3}$/dec")]):
        dd = td[td["var"] == v].reset_index(drop=True); rng = max(abs(dd["slope"])) + 1e-9
        ax.bar(range(4), dd["slope"], color=col, alpha=0.85, edgecolor="k", lw=0.6)
        for i, (sl, p) in enumerate(zip(dd["slope"], dd["p"])):
            ax.text(i, sl + np.sign(sl or 1) * rng * 0.06, star(p), ha="center",
                    va="bottom" if sl >= 0 else "top", fontsize=9, fontweight="bold")
        ax.axhline(0, color="k", lw=0.8); ax.set_xticks(range(4))
        ax.set_xticklabels([SNAME[s].replace("-", "-\n") for s in SEAS], fontsize=8)
        ax.set_title(v); ax.set_ylabel(unit); ax.set_ylim(-rng * 1.3, rng * 1.3); ax.grid(True, axis="y", ls=":", alpha=0.35)
    fig.suptitle("Seasonal trends (Sen's slope/decade;  * p<0.05, ** p<0.01, *** p<0.001, ns)",
                 fontsize=10, fontweight="bold", y=1.03)
    plt.tight_layout(); plt.savefig("Fig04_seasonal_trend.png", bbox_inches="tight"); plt.close()

    def pcorr_p(x, y, z):
        def resid(a, b):
            b = np.column_stack([np.ones_like(b), b]); return a - b @ np.linalg.lstsq(b, a, rcond=None)[0]
        r = np.corrcoef(resid(x, z), resid(y, z))[0, 1]; dof = len(x) - 3
        return r, 2 * stats.t.sf(abs(r * np.sqrt(dof / (1 - r**2))), dof)
    R = np.zeros((4, 2)); P = np.zeros((4, 2))
    for i, s in enumerate(SEAS):
        d = m[m["season"] == s].dropna(subset=["sst_C_anom", "sss_PSU_anom", "chl_log_anom"])
        x, y, c = d["sst_C_anom"].values, d["sss_PSU_anom"].values, d["chl_log_anom"].values
        R[i, 0], P[i, 0] = pcorr_p(x, c, y); R[i, 1], P[i, 1] = pcorr_p(y, c, x)
    fig, ax = plt.subplots(figsize=(5.6, 4.8))
    im = ax.imshow(R, cmap="RdBu_r", vmin=-0.7, vmax=0.7, aspect="auto")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["SST-Chl\n(control SSS)", "SSS-Chl\n(control SST)"], fontsize=9.5)
    ax.set_yticks(range(4)); ax.set_yticklabels([SNAME[s] for s in SEAS])
    for i in range(4):
        for j in range(2):
            ax.text(j, i, f"{R[i,j]:+.2f}\n{star(P[i,j])}", ha="center", va="center",
                    color="white" if abs(R[i, j]) > 0.4 else "black", fontweight="bold", fontsize=10.5)
    ax.set_title("Deseasonalized partial correlation with\nsurface Chl-a by season", pad=10)
    cb = plt.colorbar(im, ax=ax, shrink=0.85); cb.set_label("Partial correlation coefficient")
    fig.text(0.5, -0.02, "* p<0.05,  ** p<0.01,  *** p<0.001,  ns = not significant", ha="center", fontsize=8.5, style="italic")
    plt.tight_layout(); plt.savefig("Fig05_coupling_heatmap.png", bbox_inches="tight"); plt.close()

# ---------- Fig. 6: climate modes (autocorrelation-robust significance) ----------
mi = pd.read_csv("BoB_EEZ_master_with_indices.csv").sort_values("date").reset_index(drop=True)
def eff_p(a, b):
    q = pd.DataFrame({"a": a, "b": b}).dropna(); r = q["a"].corr(q["b"]); r1, r2 = q["a"].autocorr(1), q["b"].autocorr(1)
    ne = len(q) * (1 - r1 * r2) / (1 + r1 * r2); return r, 2 * stats.t.sf(abs(r * np.sqrt((ne - 2) / (1 - r * r))), ne - 2)
lags = range(8); res = {}
for nm, ic, vc in [("IOD → Chl-a", "IOD", "chl_log_anom"), ("ENSO → Chl-a", "ENSO", "chl_log_anom"), ("IOD → SSS", "IOD", "sss_PSU_anom")]:
    res[nm] = list(zip(*[eff_p(mi[ic].shift(L).values, mi[vc].values) for L in lags]))
cols = [("IOD–Chl", "IOD", "chl_log_anom"), ("ENSO–Chl", "ENSO", "chl_log_anom"), ("IOD–SSS", "IOD", "sss_PSU_anom"), ("ENSO–SSS", "ENSO", "sss_PSU_anom")]
R6 = np.zeros((4, 4)); P6 = np.zeros((4, 4))
for i, s in enumerate(SEAS):
    q = mi[mi["season"] == s]
    for j, (_, ic, vc) in enumerate(cols):
        R6[i, j] = q[[ic, vc]].dropna().corr().iloc[0, 1]
        y = q.groupby("year")[[ic, vc]].mean().dropna(); P6[i, j] = stats.pearsonr(y[ic], y[vc])[1]
fig, (a, b) = plt.subplots(1, 2, figsize=(12.5, 5), gridspec_kw={"width_ratios": [1.35, 1]})
colr = {"IOD → Chl-a": "#c8741f", "ENSO → Chl-a": "#2e8b3e", "IOD → SSS": "#127a8a"}
for nm, (rr, pp) in res.items():
    a.plot(list(lags), rr, "-", color=colr[nm], lw=2, label=nm)
    for L, r, p in zip(lags, rr, pp):
        a.plot(L, r, "o", ms=7, color=colr[nm], mfc=colr[nm] if p < 0.05 else "white", mew=1.6)
a.axhline(0, color="k", lw=0.8); a.set_xlabel("Lag (months, index leads)"); a.set_ylabel("Pearson correlation")
a.set_title("(a) Lagged correlation with climate indices", fontweight="bold", fontsize=11)
a.plot([], [], "o", color="grey", mfc="grey", label="p < 0.05 (effective n)"); a.plot([], [], "o", color="grey", mfc="white", label="not significant")
a.legend(frameon=False, fontsize=8.5, loc="upper left", bbox_to_anchor=(0.06, 0.53)); a.grid(ls=":", alpha=.5)
for s in ("top", "right"): a.spines[s].set_visible(False)
im = b.imshow(R6, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
for i in range(4):
    for j in range(4):
        st = "***" if P6[i, j] < 0.001 else "**" if P6[i, j] < 0.01 else "*" if P6[i, j] < 0.05 else ""
        b.text(j, i, minus(f"{R6[i,j]:+.2f}") + ("\n" + st if st else ""), ha="center", va="center", fontsize=9.5,
               fontweight="bold", color="white" if abs(R6[i, j]) > 0.4 else "k")
b.set_xticks(range(4)); b.set_xticklabels([c[0] for c in cols], rotation=30, ha="right")
b.set_yticks(range(4)); b.set_yticklabels(["Pre-\nmonsoon", "SW\nmonsoon", "Post-\nmonsoon", "Winter"])
b.set_title("(b) Seasonal correlation with indices", fontweight="bold", fontsize=11)
plt.colorbar(im, ax=b, fraction=0.046, pad=0.04).set_label("r")
fig.text(0.5, 0.005, "Significance: (a) filled markers, p < 0.05 with effective sample size; (b) * p < 0.05, ** p < 0.01, "
         "*** p < 0.001 from the 21 seasonal means (one value per year)", ha="center", fontsize=8.5, style="italic")
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig("Fig06_climate_modes.png", dpi=400); plt.close()

# ---------- Fig. 7: MODIS-Aqua vs Copernicus-GlobColour ----------
cv = pd.read_csv("BoB_EEZ_chl_MODIS_vs_GlobColour.csv").sort_values("date").reset_index(drop=True)
fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.4), gridspec_kw={"width_ratios": [1.75, 1]})
yrs = np.arange(2005, 2026)
for col, lab, c in [("chl_modis", "MODIS-Aqua", "#2e8b3e"), ("chl_globcolour", "Copernicus-GlobColour (merged)", "#8b3a7a")]:
    ann = cv.groupby("year")[col].mean(); r = mk.correlated_seasonal_test(cv[col].values, period=12)
    a.plot(yrs, ann.values, "-o", color=c, ms=4, lw=1.8,
           label=minus(f"{lab}: {r.slope*10:+.2f} mg m$^{{-3}}$ decade$^{{-1}}$, p = {r.p:.2f}"))
    a.plot(yrs, np.median(ann.values) + r.slope * (yrs - np.median(yrs)), "--", color=c, alpha=0.75, lw=1.3)
a.set_title("(a) Annual-mean surface Chl-a: MODIS-Aqua vs Copernicus-GlobColour", fontweight="bold", fontsize=11)
a.set_xlabel("Year"); a.set_ylabel("Annual mean Chl-a (mg m$^{-3}$)"); a.set_xticks(range(2005, 2026, 4))
a.set_xlim(2004.5, 2025.5); a.set_ylim(0.37, 0.80); a.legend(frameon=False, fontsize=9, loc="upper left"); a.grid(ls=":", alpha=.5)
q = cv.dropna(subset=["chl_modis", "chl_globcolour"]); rr = np.corrcoef(q["chl_modis"], q["chl_globcolour"])[0, 1]
b.scatter(q["chl_modis"], q["chl_globcolour"], s=12, alpha=.55, color="#4a6fa5", edgecolor="none")
lim = max(q[["chl_modis", "chl_globcolour"]].max()) * 1.05
b.plot([0, lim], [0, lim], "k--", lw=1, label="1:1 line"); b.set_xlim(0, lim); b.set_ylim(0, lim)
b.set_title(f"(b) Monthly agreement (r = {rr:.2f})", fontweight="bold", fontsize=11)
b.set_xlabel("MODIS-Aqua Chl-a (mg m$^{-3}$)"); b.set_ylabel("GlobColour Chl-a (mg m$^{-3}$)"); b.legend(frameon=False, fontsize=9); b.grid(ls=":", alpha=.5)
for ax in (a, b):
    for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("Fig07_crossvalidation.png", dpi=400); plt.close()

print("Saved Fig03-Fig07")
