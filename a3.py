import sys,pickle,json;sys.path.insert(0,'/tmp/claude-0/-home-claude/6cffcb31-050c-5a22-af60-b2eb00fe5b65/scratchpad')
import numpy as np, pandas as pd
_,_,names,w=pickle.load(open('chirps.pkl','rb')); names=list(names)
a2=json.load(open('rev/a2.json')); base4=a2['4']['base']; stab4=a2['4']['stab']
# relabel: 1 Nias, 2 Barisan/south, 3 NE coast core, 4 eastern periphery
m={2:1,1:2,3:3,4:4}; reg={n:m[base4[n]] for n in names}
mon=w.resample('MS').sum(); clim=mon.groupby(mon.index.month).mean()  # 12 x C
ann=clim.sum()
out={'region':reg,'stab':stab4}
# Oldeman
def runs(mask):
    if mask.all(): return 12
    m2=np.concatenate([mask,mask]); best=cur=0
    for v in m2:
        cur=cur+1 if v else 0; best=max(best,cur)
    return min(best,12)
old={}
for n in names:
    v=clim[n].values; wm=runs(v>200); dm=runs(v<100)
    T='A' if wm>9 else 'B' if wm>=7 else 'C' if wm>=5 else 'D' if wm>=3 else 'E'
    s=1 if dm<2 else 2 if dm<=3 else 3 if dm<=6 else 4
    old[n]=dict(wet=int(wm),dry=int(dm),type=f'{T}{s}')
out['oldeman']=old
# BMKG dasarian criteria on climatological dasarian means
d=w.copy(); d['y']=d.index.year; d['m']=d.index.month; d['das']=np.minimum((d.index.day-1)//10,2)
das=d.groupby(['y','m','das'])[names].sum().groupby(level=[1,2]).mean()  # 36 x C
bm={}
for n in names:
    v=das[n].values; below=v<50
    # longest circular run of dasarian <50
    m2=np.concatenate([below,below]); best=cur=0
    for x in m2: cur=cur+1 if x else 0; best=max(best,cur)
    bm[n]=dict(n_below50=int(below.sum()),longest_dry_run=int(min(best,36)),min_das=round(float(v.min()),1))
out['bmkg']=bm
# region summaries
rows=[]
for r in (1,2,3,4):
    mem=[n for n in names if reg[n]==r]; prof=clim[mem].mean(1)
    rows.append(dict(region=r,n=len(mem),members=mem,annual=round(float(ann[mem].mean())),ann_min=(ann[mem].idxmin(),round(float(ann[mem].min()))),ann_max=(ann[mem].idxmax(),round(float(ann[mem].max()))),
                     monthly=[round(float(x)) for x in prof.values],oldeman=sorted(set(old[x]['type'] for x in mem)),
                     bmkg_dry_units=sum(bm[x]['longest_dry_run']>=3 for x in mem),stable=[x for x in mem if stab4[x]>=0.75]))
out['regions']=rows
for r in rows: print(r)
for n in names: print(n,reg[n],stab4[n],old[n],bm[n],round(float(ann[n])))
# pixel counts from boundary geojson
g=json.load(open('/mnt/user-data/attachments/6cffcb31-050c-5a22-af60-b2eb00fe5b65/6ec60b52-zona_curah_hujan_sumut.geojson'))
def polys(geom):
    t=geom['type']
    if t=='Polygon': yield geom['coordinates']
    elif t=='MultiPolygon':
        for p in geom['coordinates']: yield p
    elif t=='GeometryCollection':
        for gg in geom['geometries']: yield from polys(gg)
def area(ring):
    x=np.array([p[0] for p in ring]);y=np.array([p[1] for p in ring]); return 0.5*abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))
REN={'Labuhan Batu':'Labuhanbatu','Labuhan Batu Selatan':'Labuhanbatu Selatan','Labuhan Batu Utara':'Labuhanbatu Utara','Toba Samosir':'Toba','Kota Tanjung Balai':'Kota Tanjungbalai'}
pix={}
for f in g['features']:
    nm=REN.get(f['properties']['nama'],f['properties']['nama']); a=0
    for p in polys(f['geometry']):
        a+=area(p[0])-sum(area(h) for h in p[1:])
    pix[nm]=round(a/0.05**2,1)
out['pixels']=pix
print(sorted(pix.items(),key=lambda x:x[1])[:10])
json.dump(out,open('rev/a3.json','w'),indent=1,ensure_ascii=False)
