# Changelog

## Unreleased

### Changed
- `03_trend_analysis.py`: the continuous December–February winter sensitivity
  test (section C), previously only printed to the screen, is now also saved to
  `tables/TableS1_winter_DecFeb_sensitivity.csv` (Supplementary Table S1, rows
  marked †; manuscript §3.2). The computation is unchanged, and
  `Table2_annual_trends.csv` and `TableS1_seasonal_trends.csv` are reproduced
  identically.
- `README.md` and `tables/README.md`: new output file listed.

## 1.2.0 — revisions in response to review

### Added
- `03b_seasonal_autocorrelation.py` and `tables/TableS1b_seasonal_autocorrelation.csv`:
  lag-1 autocorrelation of each detrended seasonal series used in the
  per-season Mann–Kendall tests (three variables × four seasons, 21 yearly
  values each). All values satisfy |r1| ≤ 0.25, below the approximate 95% limit
  of ±0.43, so no correction for inter-annual autocorrelation is applied
  (manuscript §2.5).

### Changed
- `04_coupling_analysis.py`: the coverage-threshold test now also reports
  thresholds of 10% and 20% alongside the main 15% threshold;
  `tables/TableS2_seasonal_coupling.csv` gains the corresponding columns
  (manuscript §2.6, §3.4). The existing columns (full record and 15%) are
  unchanged, and Tables S3 and S4a are reproduced identically.
- README, `tables/README.md` and `code/README.md`: run order, table list and
  result descriptions updated; terminology aligned with the revised manuscript
  ("dominant partial association"; significance "after false-discovery-rate
  correction").
- `CITATION.cff`: version 1.2.0.

## 1.1.0 — aligned with the revised manuscript

### Fixed
- **Trend test (scripts 03, 06).** Missing Chl-a months were dropped before the
  seasonal Mann–Kendall test, which shifted every later month into the wrong
  calendar position. Missing months are now kept as NaN.
- **Serial dependence (03, 06).** Full-record trends now use the Hirsch–Slack
  form of the seasonal Mann–Kendall test, as reported in the paper. The
  uncorrected p-value is still saved for comparison.
- **HTTPS verification (05).** Certificate checking was switched off when
  downloading the NOAA indices; it is now on.
- **Earth Engine asset (00).** The asset path is now a clearly marked
  `EEZ_ASSET` setting (default: the authors' asset), with instructions for
  other users.

### Added
- `00b_get_EEZ.py`: downloads the EEZ (Marine Regions v12, MRGID 8481) and the
  Bay of Bengal polygon, and writes `EEZ.shp` and `EEZ.zip` for Earth Engine.
- Seasonal-mean robustness test (04, Table S4a) and offshore sub-domain test
  (`04b_offshore_test.py`, Table S4b).
- Effective-sample-size and seasonal-mean significance for the climate-mode
  correlations (05, `climate_mode_correlations.csv`).
- Continuous December–February winter sensitivity test (03).
- False-discovery-rate field significance for the trend maps (08).
- Code for Figs. 1–2 (09), and Figs. 6–7 (07), which previously had no script.
- `10_supporting_numbers.py`: reproduces the remaining numbers quoted in the
  text (estuarine-cell share, latitude-weighting effect, GLORYS SST trend,
  record statistics, Chl-a dipole by zone).
- `Table2_annual_trends.csv`, `CHANGELOG.md`, pinned package versions.
- Author, affiliation and contact details in README, `CITATION.cff` and `LICENSE`.

### Changed
- "OC-CCI" renamed to Copernicus-GlobColour throughout (the product actually
  used is CMEMS OCEANCOLOUR_GLO_BGC_L4_MY_009_104); script 06 and its files
  renamed accordingly.
- Figures renamed to match the paper's numbering (Fig01–Fig11).
- README: record length stated as 21 years; results, data table and run order
  updated; the "where each result comes from" table added.

## 1.0.0 — initial version
