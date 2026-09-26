/* ============================================================================
   PHASE 1 — Bangladesh EEZ monthly time series (Google Earth Engine, JS)
   ----------------------------------------------------------------------------
   Variables in THIS script:
     • SST  : NOAA OISST v2.1        (NOAA/CDR/OISST/V2_1),  0.25°, daily
     • Chl-a: MODIS-Aqua L3SMI       (NASA/OCEANDATA/MODIS-Aqua/L3SMI), 4 km, daily
   Output:
     • One CSV to Google Drive: monthly EEZ-averaged SST & Chl-a (2005–2025)
       with valid-pixel coverage + season tag.
   NOTE:
     • SSS (GLORYS12V1) is NOT here — GEE lacks a good coastal SSS product.
       That comes from Copernicus in a separate script (Phase 1b).
     • MODIS-Aqua chl starts Jul-2002; period here = 2005–2025 (21 full years).
   ========================================================================== */


/* ----------------------------------------------------------------------------
   1. STUDY AREA — Bangladesh EEZ  (using your uploaded asset)
--------------------------------------------------------------------------- */
// Bangladesh EEZ polygon (Marine Regions, Maritime Boundaries v12, MRGID 8481),
// uploaded as an Earth Engine table asset (EEZ.zip from code/00b_get_EEZ.py).
// The path below is the authors' asset. Other users: upload EEZ.zip to your own
// Earth Engine project and replace this path with your asset ID.
var EEZ_ASSET = 'projects/xenon-depth-502316-d8/assets/BD_EEZ';
var eez = ee.FeatureCollection(EEZ_ASSET);

/* Test box kept below only for reference — leave commented out.
var eezApprox = ee.Geometry.Polygon(
  [[[88.0, 22.6], [92.5, 22.6], [92.5, 15.5], [88.0, 15.5], [88.0, 22.6]]], null, false);
var eez = ee.FeatureCollection([ee.Feature(eezApprox)]);
*/

var roi = eez.geometry();
Map.centerObject(roi, 6);
Map.addLayer(roi, {color: 'red'}, 'Bangladesh EEZ');


/* ----------------------------------------------------------------------------
   2. PARAMETERS
--------------------------------------------------------------------------- */
var startYear = 2005;   // start of analysis period (MODIS-Aqua chl begins Jul-2002)
var endYear   = 2025;   // last complete year (2005–2025 inclusive = 21 full years)
//                         → for exactly 20 years, set endYear = 2024.
var startDate = ee.Date.fromYMD(startYear, 1, 1);
var endDate   = ee.Date.fromYMD(endYear + 1, 1, 1);   // exclusive upper bound

var SST_SCALE = 25000;  // ~0.25° in meters, for reduceRegion on OISST
var CHL_SCALE = 4616;   // MODIS L3SMI native pixel size (m)


/* ----------------------------------------------------------------------------
   3. SOURCE COLLECTIONS
--------------------------------------------------------------------------- */
// OISST sst band is stored scaled; multiply by 0.01 to get °C.
var oisst = ee.ImageCollection('NOAA/CDR/OISST/V2_1')
              .filterDate(startDate, endDate)
              .select('sst');

// MODIS-Aqua chlorophyll-a (mg/m3); scale/offset already applied in GEE.
var modis = ee.ImageCollection('NASA/OCEANDATA/MODIS-Aqua/L3SMI')
              .filterDate(startDate, endDate)
              .select('chlor_a');


/* ----------------------------------------------------------------------------
   4. REFERENCE PIXEL COUNTS (to compute % valid coverage per month)
      Max possible pixels inside the EEZ at each product's resolution.
--------------------------------------------------------------------------- */
var maxChlPix = ee.Number(
  ee.Image.constant(1).reduceRegion({
    reducer: ee.Reducer.count(), geometry: roi,
    scale: CHL_SCALE, maxPixels: 1e13, bestEffort: true
  }).get('constant'));


/* ----------------------------------------------------------------------------
   5. SEASON LOOKUP (BoB monsoon-based, not generic DJF/MAM)
      pre_monsoon = Mar–May | sw_monsoon = Jun–Sep
      post_monsoon = Oct–Nov | winter (NE monsoon) = Dec–Feb
--------------------------------------------------------------------------- */
var seasonDict = ee.Dictionary({
  '1':'winter','2':'winter','3':'pre_monsoon','4':'pre_monsoon','5':'pre_monsoon',
  '6':'sw_monsoon','7':'sw_monsoon','8':'sw_monsoon','9':'sw_monsoon',
  '10':'post_monsoon','11':'post_monsoon','12':'winter'});


/* ----------------------------------------------------------------------------
   6. MONTHLY LOOP  → one feature per month
--------------------------------------------------------------------------- */
var nMonths  = endDate.difference(startDate, 'month').round();
var monthSeq = ee.List.sequence(0, nMonths.subtract(1));

var MISSING = -999;   // sentinel for fully-masked months → treat as NaN in analysis
//  (EE drops null-valued keys from a Dictionary, so a real number is used
//   as the default; this guarantees .get() below never hits a missing key.)

var monthly = ee.FeatureCollection(monthSeq.map(function (m) {
  m = ee.Number(m);
  var mStart = startDate.advance(m, 'month');
  var mEnd   = mStart.advance(1, 'month');
  var year   = mStart.get('year');
  var month  = mStart.get('month');

  // ---- SST: monthly mean composite → area-mean over EEZ (°C) ----
  var sstImg  = oisst.filterDate(mStart, mEnd).mean().multiply(0.01);
  var sstDict = ee.Dictionary({sst: MISSING}).combine(
    sstImg.reduceRegion({
      reducer: ee.Reducer.mean(), geometry: roi,
      scale: SST_SCALE, maxPixels: 1e13, bestEffort: true
    }), true);   // true = real value overwrites the sentinel when present

  // ---- Chl-a: monthly composite (mg/m3) ----
  var chlImg = modis.filterDate(mStart, mEnd).mean();

  // (a) mean of log10(chl); geometric mean = 10^chl_logmean (done in analysis)
  var logDict = ee.Dictionary({chlor_a: MISSING}).combine(
    chlImg.log10().reduceRegion({
      reducer: ee.Reducer.mean(), geometry: roi,
      scale: CHL_SCALE, maxPixels: 1e13, bestEffort: true
    }), true);

  // (b) spatial median over EEZ (mg/m3)
  var medDict = ee.Dictionary({chlor_a: MISSING}).combine(
    chlImg.reduceRegion({
      reducer: ee.Reducer.median(), geometry: roi,
      scale: CHL_SCALE, maxPixels: 1e13, bestEffort: true
    }), true);

  // (c) valid pixel count + % coverage (KEY for monsoon cloud gaps)
  var cntDict = ee.Dictionary({chlor_a: 0}).combine(
    chlImg.reduceRegion({
      reducer: ee.Reducer.count(), geometry: roi,
      scale: CHL_SCALE, maxPixels: 1e13, bestEffort: true
    }), true);
  var chlCnt = ee.Number(cntDict.get('chlor_a'));
  var chlPct = chlCnt.divide(maxChlPix).multiply(100);

  var season = seasonDict.get(ee.Number(month).format('%d'));

  return ee.Feature(null, {
    date:          mStart.format('YYYY-MM'),
    year:          year,
    month:         month,
    season:        season,
    sst_C:         sstDict.get('sst'),
    chl_logmean:   logDict.get('chlor_a'),   // geometric mean = 10^chl_logmean
    chl_median:    medDict.get('chlor_a'),
    chl_validpix:  chlCnt,
    chl_validpct:  chlPct
  });
}));


/* ----------------------------------------------------------------------------
   7. QUICK LOOK (optional) — print a chart in the console
--------------------------------------------------------------------------- */
print('First 12 monthly records:', monthly.limit(12));

print(ui.Chart.feature.byFeature(monthly, 'date', ['sst_C'])
        .setChartType('LineChart')
        .setOptions({title: 'Monthly EEZ-mean SST (°C)', legend: {position:'none'}}));

print(ui.Chart.feature.byFeature(monthly, 'date', ['chl_median'])
        .setChartType('LineChart')
        .setOptions({title: 'Monthly EEZ median Chl-a (mg/m3)', legend:{position:'none'}}));

print(ui.Chart.feature.byFeature(monthly, 'date', ['chl_validpct'])
        .setChartType('ColumnChart')
        .setOptions({title: 'Chl-a valid-pixel coverage (%) — watch monsoon months',
                     legend: {position:'none'}}));


/* ----------------------------------------------------------------------------
   8. EXPORT → Google Drive (CSV)
--------------------------------------------------------------------------- */
Export.table.toDrive({
  collection: monthly,
  description: 'BoB_EEZ_SST_Chl_monthly_2005_2025',
  fileNamePrefix: 'BoB_EEZ_SST_Chl_monthly_2005_2025',
  fileFormat: 'CSV',
  selectors: ['date','year','month','season',
              'sst_C','chl_logmean','chl_median','chl_validpix','chl_validpct']
});

/* After running: open the 'Tasks' tab (top-right) → click RUN on the export.
   The CSV lands in your Google Drive root. */
