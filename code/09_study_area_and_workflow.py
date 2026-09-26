"""
============================================================================
09_study_area_and_workflow.py
----------------------------------------------------------------------------
Purpose : Fig01_study_area.png (study-area map) and Fig02_workflow.png
          (methodological workflow diagram).
Input   : EEZ.shp and bob_iho.geojson (from 00b_get_EEZ.py)
Note    : the first run of cartopy downloads Natural Earth coastlines.
============================================================================
"""
import numpy as np
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
from matplotlib.patches import Rectangle, FancyBboxPatch, FancyArrowPatch

# ---------------- Fig. 1: study area ----------------
eez = gpd.read_file('EEZ.shp'); bob = gpd.read_file('bob_iho.geojson')
fig=plt.figure(figsize=(7.2,8)); ax=plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([78,96,5,25],ccrs.PlateCarree())
ax.add_geometries(bob.geometry,ccrs.PlateCarree(),facecolor='#eaf2f8',edgecolor='none',zorder=0)
ax.add_feature(cf.LAND.with_scale('10m'),facecolor='#ece8e0',edgecolor='none',zorder=1)
ax.add_feature(cf.COASTLINE.with_scale('10m'),lw=0.5,edgecolor='#6b7b86',zorder=2)
ax.add_feature(cf.BORDERS.with_scale('10m'),lw=0.4,edgecolor='#9a9a9a',linestyle=':',zorder=2)
ax.add_geometries(eez.geometry,ccrs.PlateCarree(),facecolor='#c0392b',alpha=0.35,edgecolor='#7b1e2b',lw=1.2,zorder=3)
ax.annotate('Bangladesh\nEEZ',xy=(90.6,19.6),xytext=(93.3,17.3),fontsize=10,fontweight='bold',color='#7b1e2b',ha='center',
            arrowprops=dict(arrowstyle='-',color='#7b1e2b',lw=0.8),zorder=5)
ax.text(86.2,12.2,'Bay of Bengal',fontsize=12,style='italic',color='#2e5f7a',ha='center')
for x,y,t in [(79.8,21.5,'INDIA'),(94.2,21.8,'MYANMAR'),(90.2,23.9,'BANGLADESH')]:
    ax.text(x,y,t,fontsize=8,color='#555',ha='center',zorder=4)
gl=ax.gridlines(draw_labels=True,lw=0.3,color='grey',alpha=0.5,linestyle=':'); gl.top_labels=gl.right_labels=False
# scale bar (200 km at ~7N)
lat0=7.0; x0=86.3; km=200; dlon=km/(111.32*np.cos(np.deg2rad(lat0)))
ax.add_patch(Rectangle((x0,lat0),dlon,0.18,facecolor='k',transform=ccrs.PlateCarree(),zorder=6))
ax.add_patch(Rectangle((x0+dlon,lat0),dlon,0.18,facecolor='white',edgecolor='k',lw=0.6,transform=ccrs.PlateCarree(),zorder=6))
for k,xx in enumerate([x0,x0+dlon,x0+2*dlon]): ax.text(xx,lat0+0.35,f'{k*km}',fontsize=7,ha='center',zorder=6)
ax.text(x0+2*dlon+0.3,lat0+0.05,'km',fontsize=7,zorder=6)
ax.annotate('N',xy=(95.2,24.4),xytext=(95.2,22.9),ha='center',fontsize=10,fontweight='bold',arrowprops=dict(arrowstyle='-|>',color='k',lw=1.2),zorder=6)
ax.set_title('Study area: Bangladesh EEZ, northern Bay of Bengal',fontweight='bold',fontsize=11)
plt.savefig('Fig01_study_area.png', dpi=400, bbox_inches='tight'); plt.close()

# ---------------- Fig. 2: workflow ----------------
fig,ax=plt.subplots(figsize=(12,10.6)); ax.set_xlim(0,12); ax.set_ylim(0,10.6); ax.axis('off')
C={'acq':('#2c6e9b','#e8f1f8'),'pre':('#3a7d3a','#edf5ea'),'ana':('#c47a2c','#fdf1e6'),'out':('#8a3f7d','#f5e9f3')}
def box(x,y,w,h,txt,st,bold=False,fs=9.2):
    ec,fc=C[st]; ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.02,rounding_size=0.08',ec=ec,fc=fc,lw=1.8))
    ax.text(x+w/2,y+h/2,txt,ha='center',va='center',fontsize=fs,fontweight='bold' if bold else 'normal',linespacing=1.3)
def frame(y,h,title,st):
    ec,_=C[st]; ax.add_patch(FancyBboxPatch((0.15,y),11.7,h,boxstyle='round,pad=0.02,rounding_size=0.05',ec=ec,fc='none',lw=1,ls='--',alpha=0.6))
    ax.text(0.3,y+h+0.1,title,ha='left',fontsize=12.5,fontweight='bold',color=ec)
def arr(x1,y1,x2,y2,ls='-'):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=13,lw=1.4,color='#333',ls=ls))
# 1 acquisition
frame(8.55,1.45,'DATA ACQUISITION','acq')
acq=[('SST\nNOAA OISST v2.1\n0.25°, daily (GEE)'),('Chl-a\nMODIS-Aqua L3SMI\n4 km, daily (GEE)'),('SSS (and SST for maps)\nGLORYS12V1\n1/12°, monthly (CMEMS)'),
     ('Chl-a (validation, maps)\nCopernicus-GlobColour L4\n4 km, monthly (CMEMS)'),('Climate indices\nDMI and Niño 3.4\nmonthly (NOAA PSL)')]
w=2.18
for k,t in enumerate(acq): box(0.3+k*(w+0.12),8.68,w,1.2,t,'acq',fs=8.6)
# 2 pre-processing
frame(6.25,1.7,'PRE-PROCESSING','pre')
box(1.6,7.05,8.8,0.8,'Mask to the Bangladesh EEZ polygon (Marine Regions)  →  cos-latitude-weighted monthly EEZ mean\nChl-a: geometric (log$_{10}$) mean and valid-pixel percentage recorded','pre',fs=9.2)
box(2.6,6.38,6.8,0.5,'Master monthly dataset: 252 months (2005–2025) of EEZ-mean SST, SSS and Chl-a','pre',bold=True,fs=9.4)
for k in range(5): arr(0.3+k*(w+0.12)+w/2,8.68,6,7.87)
arr(6,7.05,6,6.9)
# 3 analysis
frame(2.75,2.95,'ANALYSIS','ana')
box(0.6,4.75,5.1,0.75,'Deseasonalisation\nanomaly = value − monthly climatology (2005–2025)','ana')
box(6.3,4.75,5.1,0.75,'Seasonal split\npre-monsoon · SW monsoon · post-monsoon · winter','ana')
arr(6,6.38,3.15,5.52); arr(6,6.38,8.85,5.52)
an=['Trends\nHirsch–Slack seasonal MK\n+ seasonal Sen’s slope;\nMK on seasonal means;\nper-pixel maps',
    'Coupling\nseasonal partial\ncorrelation\n(SST | SSS, SSS | SST)',
    'Climate modes\nlagged and seasonal\ncorrelation with\nDMI and Niño 3.4',
    'Robustness\ncoverage threshold (15%);\noffshore sub-domain;\nseasonal-mean test;\nGlobColour cross-check']
w2=2.72
for k,t in enumerate(an): box(0.35+k*(w2+0.13),2.9,w2,1.55,t,'ana',fs=8.8)
xs=[0.35+k*(w2+0.13)+w2/2 for k in range(4)]
arr(3.15,4.75,xs[0],4.47); arr(3.15,4.75,xs[1],4.47,'--'); arr(3.15,4.75,xs[2],4.47,'--')
arr(8.85,4.75,xs[1],4.47,'--'); arr(8.85,4.75,xs[2],4.47,'--'); arr(8.85,4.75,xs[3],4.47)
# 4 outputs
frame(0.2,1.9,'OUTPUTS','out')
outs=['Trends\nFigs. 3–4, Table 2','Spatial patterns\nFigs. 8–11','Driver handover\nFig. 5, Table 3','Climate modes\nFig. 6','Robustness and\nvalidation\nFig. 7, Suppl. Tables']
for k,t in enumerate(outs): box(0.3+k*(w+0.12),0.35,w,1.25,t,'out',fs=9)
ox=[0.3+k*(w+0.12)+w/2 for k in range(5)]
arr(xs[0],2.9,ox[0],1.62); arr(xs[0],2.9,ox[1],1.62); arr(xs[1],2.9,ox[2],1.62); arr(xs[2],2.9,ox[3],1.62); arr(xs[3],2.9,ox[4],1.62)
plt.savefig('Fig02_workflow.png', dpi=400, bbox_inches='tight', facecolor='white'); plt.close()

print("Saved Fig01_study_area.png and Fig02_workflow.png")
