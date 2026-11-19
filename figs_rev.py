"""Revised figures (k = 4 primary scheme) for the ICDSOS revision."""
import json, pickle, sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.patheffects as pe
from matplotlib.patches import Polygon as MPoly, Patch
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
S='/tmp/claude-0/-home-claude/6cffcb31-050c-5a22-af60-b2eb00fe5b65/scratchpad/'
OUT='/home/claude/hasil_revisi/'
a1=json.load(open(S+'rev/a1.json')); a2=json.load(open(S+'rev/a2.json')); a3=json.load(open(S+'rev/a3.json'))
A=pickle.load(open(S+'rev/a1.pkl','rb')); _,_,names,w=pickle.load(open(S+'chirps.pkl','rb')); names=list(names)
COL={1:"#2a78d6",2:"#eb6834",3:"#1baf7a",4:"#eda100"}; MK={1:"o",2:"s",3:"^",4:"D"}
INK,MUTED,GRID="#1f1f1e","#6b6a63","#e4e3de"; HALO=[pe.withStroke(linewidth=2.5,foreground="white")]
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.edgecolor":MUTED,"axes.labelcolor":INK,"xtick.color":MUTED,"ytick.color":MUTED,
  "axes.spines.top":False,"axes.spines.right":False,"axes.grid":True,"grid.color":GRID,"grid.linewidth":0.6,"savefig.dpi":300,"savefig.bbox":"tight"})
reg={n:int(r) for n,r in a3['region'].items()}; stab=a2['4']['stab']
RN={1:"Nias archipelago",2:"Barisan highlands and south",3:"North-eastern coastal plain",4:"Eastern fringe"}
# Fig 3 validation
tab=a1['table']; ser=[('MSSA','MSSA seasonal (RC1–4)'),('Raw','Raw daily'),('Clim','Monthly climatology'),('Harmonic2','Harmonic regression (2 harmonics)')]
fig,axes=plt.subplots(1,2,figsize=(9.5,3.8))
for met,ax,lab in [('sil',axes[0],'Average silhouette'),('ch',axes[1],'Calinski–Harabasz index')]:
    for i,(key,lg) in enumerate(ser):
        d=[r for r in tab if r['input']==key]
        ax.plot([r['k'] for r in d],[r[met] for r in d],color=COL[i+1],lw=2,marker=MK[i+1],ms=6,markeredgecolor="white",markeredgewidth=1.3,label=lg)
    ax.axvline(4,color=MUTED,lw=1,ls=":"); ax.set_xticks(range(2,9)); ax.set_xlabel("Number of clusters (k)"); ax.set_ylabel(lab)
    ax.set_title(lab,fontsize=10,loc="left")
axes[0].legend(frameon=False,fontsize=8,loc="lower left",bbox_to_anchor=(0,1.08),ncol=2)
fig.tight_layout(); fig.savefig(OUT+'Fig3_cluster_validation.png'); plt.close(fig)
# Fig 4 dendrogram
F=A['main']; F=(F-F.mean(1,keepdims=True))/F.std(1,ddof=1,keepdims=True); Lk=linkage(F,'ward'); C=len(names)
leafreg={i:reg[n] for i,n in enumerate(names)}; node=dict(leafreg)
for i,(a,b,_,_) in enumerate(Lk):
    za,zb=node[int(a)],node[int(b)]; node[C+i]=za if za==zb else None
fig,ax=plt.subplots(figsize=(9.5,4.4))
dendrogram(Lk,labels=[n for n in names],leaf_rotation=90,leaf_font_size=7.5,link_color_func=lambda k:COL[node[k]] if node.get(k) else "#9a998f",ax=ax)
for t in ax.get_xticklabels():
    t.set_color(INK); n=t.get_text()
    if stab[n]<0.75: t.set_text(n+" *"); t.set_fontstyle('italic')
ax.set_xticklabels([t.get_text() for t in ax.get_xticklabels()])
cut=(Lk[-3,2]+Lk[-4,2])/2; ax.axhline(cut,color=MUTED,lw=1,ls=":"); ax.grid(False); ax.set_ylabel("Ward linkage distance")
ax.set_title("Ward dendrogram of MSSA seasonal signals, cut at k = 4 (* transitional unit)",fontsize=10,loc="left")
fig.savefig(OUT+'Fig4_dendrogram.png'); plt.close(fig)
# Fig 5 seasonal profiles
months=["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
fig,ax=plt.subplots(figsize=(7.8,4.2))
for r in a3['regions']:
    z=r['region']; ax.plot(range(1,13),r['monthly'],color=COL[z],lw=2,marker=MK[z],ms=6,markeredgecolor="white",markeredgewidth=1.3,label=f"Region {z} – {RN[z]} (n = {r['n']})")
ax.axhline(100,color=MUTED,lw=0.8,ls="--"); ax.axhline(200,color=MUTED,lw=0.8,ls="--")
ax.text(12.35,100,"Oldeman dry < 100 mm",fontsize=7.5,color=MUTED,va="center"); ax.text(12.35,200,"Oldeman wet > 200 mm",fontsize=7.5,color=MUTED,va="center")
ax.set_xticks(range(1,13)); ax.set_xticklabels(months); ax.set_xlim(0.6,12.3); ax.set_ylim(60,450); ax.set_axisbelow(True)
ax.set_ylabel("Mean monthly precipitation (mm)"); ax.legend(frameon=False,fontsize=8,loc="upper left")
ax.set_title("Seasonal rainfall regime by region (CHIRPS 1995–2024)",fontsize=10,loc="left")
fig.savefig(OUT+'Fig5_seasonal_profile.png'); plt.close(fig)
# Fig 6 map
g=json.load(open('/mnt/user-data/attachments/6cffcb31-050c-5a22-af60-b2eb00fe5b65/6ec60b52-zona_curah_hujan_sumut.geojson'))
REN={'Labuhan Batu':'Labuhanbatu','Labuhan Batu Selatan':'Labuhanbatu Selatan','Labuhan Batu Utara':'Labuhanbatu Utara','Toba Samosir':'Toba','Kota Tanjung Balai':'Kota Tanjungbalai'}
def polys(geom):
    t=geom['type']
    if t=='Polygon': yield geom['coordinates']
    elif t=='MultiPolygon':
        for p in geom['coordinates']: yield p
    elif t=='GeometryCollection':
        for gg in geom['geometries']: yield from polys(gg)
plt.rcParams.update({"font.size":9,"axes.grid":False})
fig,ax=plt.subplots(figsize=(7.2,7.2)); cent={}
for f in g['features']:
    nm=REN.get(f['properties']['nama'],f['properties']['nama']); z=reg[nm]; pts=[]
    for p in polys(f['geometry']):
        ring=np.array(p[0]); pts.append(ring)
        ax.add_patch(MPoly(ring,closed=True,facecolor=COL[z],edgecolor="white",linewidth=0.6))
        if stab[nm]<0.75: ax.add_patch(MPoly(ring,closed=True,facecolor="none",edgecolor="#ffffff",hatch="////",linewidth=0))
    cent[nm]=max(pts,key=len).mean(0)
ax.set_aspect(1/np.cos(np.deg2rad(2.3))); ax.set_xlim(96.85,100.65); ax.set_ylim(-0.65,4.45)
ax.set_xticks(range(97,101)); ax.set_yticks(range(0,5))
ax.set_xticklabels([f"{x}°E" for x in range(97,101)],color=MUTED); ax.set_yticklabels([f"{y}°N" for y in range(0,5)],color=MUTED)
for s in ax.spines.values(): s.set_visible(True); s.set_color("#b9b8b0"); s.set_linewidth(0.6)
ax.set_facecolor("#f4f7fa")
for z,(x,y) in {1:(97.55,0.85),2:(99.05,1.15),3:(99.05,3.42),4:(100.15,2.2)}.items():
    ax.text(x,y,str(z),fontsize=16,fontweight="bold",color="white",ha="center",va="center",path_effects=[pe.withStroke(linewidth=3,foreground=INK)])
for nm,lbl,dx,dy in [("Kota Medan","Medan",0.25,0.18),("Kota Sibolga","Sibolga",-0.55,-0.05),("Kota Gunungsitoli","Gunungsitoli",-0.15,0.32),("Kota Padangsidimpuan","Padangsidimpuan",0.75,0.05)]:
    x,y=cent[nm]; ax.annotate(lbl,(x,y),(x+dx,y+dy),fontsize=7.5,color=INK,ha="center",path_effects=HALO,arrowprops=dict(arrowstyle="-",color=MUTED,lw=0.5)); ax.plot(x,y,"o",ms=2.5,color=INK)
ax.text(99.05,2.42,"Lake Toba",fontsize=6.5,color=MUTED,ha="left",va="center",style="italic",path_effects=HALO)
ax.text(97.4,0.25,"Indian Ocean",fontsize=8,color=MUTED,style="italic",ha="center")
ax.text(99.75,3.75,"Strait of Malacca",fontsize=8,color=MUTED,style="italic",ha="center",rotation=-38)
km=1/(111.32*np.cos(np.deg2rad(0.5))); x0,y0=99.55,-0.45
ax.plot([x0,x0+50*km],[y0,y0],color=INK,lw=2,solid_capstyle="butt"); ax.plot([x0+50*km,x0+100*km],[y0,y0],color=MUTED,lw=2,solid_capstyle="butt")
for d,l in ((0,"0"),(50,"50"),(100,"100 km")): ax.text(x0+d*km,y0+0.07,l,fontsize=7,ha="center",color=INK)
ax.annotate("",xy=(100.45,4.3),xytext=(100.45,3.95),arrowprops=dict(arrowstyle="-|>",color=INK,lw=1.2)); ax.text(100.45,4.33,"N",ha="center",va="bottom",fontsize=9,color=INK,fontweight="bold")
nreg={r['region']:r['n'] for r in a3['regions']}
ax.legend(handles=[Patch(facecolor=COL[z],edgecolor="white",label=f"Region {z} – {RN[z]} (n = {nreg[z]})") for z in (1,2,3,4)]+
          [Patch(facecolor="#9a998f",edgecolor="white",hatch="////",label="Transitional unit (kept its region in < 75% of the 15 alternative configurations)")],
          loc="upper center",bbox_to_anchor=(0.5,-0.06),frameon=False,fontsize=8)
ax.set_title("Rainfall regions of North Sumatra (MSSA RC1–4 + Ward, k = 4)",fontsize=10,loc="left",color=INK)
fig.text(0.13,-0.07,"Data: CHIRPS daily 1995–2024; boundaries: geoBoundaries ADM2 and FAO GAUL 2015 (Nias Selatan).",fontsize=6.5,color=MUTED)
fig.savefig(OUT+'Fig6_zone_map.png'); plt.close(fig)
print('ok')
