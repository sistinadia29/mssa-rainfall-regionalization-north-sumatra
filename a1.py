import sys,pickle,json,time;sys.path.insert(0,'/tmp/claude-0/-home-claude/6cffcb31-050c-5a22-af60-b2eb00fe5b65/scratchpad')
from rssa_like import *
from sklearn.metrics import adjusted_rand_score as ari
OUT={}
_,_,names,w=pickle.load(open('chirps.pkl','rb')); names=list(names)
P=w.values;N,C=P.shape
def stdz(X): return (X-X.mean(0))/X.std(0,ddof=1)
def rowz(F): return (F-F.mean(1,keepdims=True))/F.std(1,ddof=1,keepdims=True)
def per(x): x=x-x.mean();f=np.abs(np.fft.rfft(x))[1:];return len(x)/(f.argmax()+1)
def wcor_fn(R,L):
    n=R[0].shape[0];K=n-L+1;wt=np.array([min(t+1,L,K,n-t) for t in range(n)],float);ww=np.tile(wt,R[0].shape[1])
    fl=[r.ravel('F') for r in R]
    return lambda i,j:(fl[i]*fl[j]*ww).sum()/np.sqrt((fl[i]**2*ww).sum()*(fl[j]**2*ww).sum())
def ward(F,k): return fcluster(linkage(rowz(F),'ward'),k,'maxclust')
def score(F,k): l=ward(F,k); return l, silhouette_score(rowz(F),l), calinski_harabasz_score(rowz(F),l)
Z=stdz(P)
t0=time.time(); sh,R=mssa(Z,365,10); print('mssa',round(time.time()-t0))
main=sum(R[:4]).T; base={k:ward(main,k) for k in range(2,9)}
pickle.dump((sh,R),open('rev/mssa365.pkl','wb'))
# component selection per config
comp={}
def describe(tag,shr,Rr,L,nk=8):
    ps=[per(r.mean(1)) for r in Rr[:nk]]; wc=wcor_fn(Rr[:nk],L)
    sel=[j for j,p in enumerate(ps) if 160<p<200 or 330<p<400]
    pairs=[]
    for j in range(nk-1):
        pairs.append(round(abs(wc(j,j+1)),2))
    comp[tag]=dict(periods=[round(p,1) for p in ps],share=[round(100*s,2) for s in shr[:nk]],selected=[j+1 for j in sel],
                   adjacent_wcor=pairs, max_cross=round(max([abs(wc(i,j)) for i in sel for j in range(nk) if j not in sel] or [0]),2))
    print(tag,comp[tag])
    return sel
describe('L=365',sh,R,365,10)
for L in [180,730,1095]:
    s2,R2=mssa(Z,L,8); sel=describe(f'L={L}',s2,R2,L)
yrs=w.index.year.values
for nm,msk in [('1995-2009',yrs<=2009),('2010-2024',yrs>=2010)]:
    s2,R2=mssa(stdz(P[msk]),365,8); describe(nm,s2,R2,365)
OUT['components']=comp
# 122-day inclusion
ps=[per(r.mean(1)) for r in R]; idx122=[j for j,p in enumerate(ps) if 110<p<135]; print('122d comps',[j+1 for j in idx122])
F6=sum(R[j] for j in [0,1,2,3]+idx122).T
OUT['with122']={k:dict(sil=round(score(F6,k)[1],3),ari=round(ari(base[k],ward(F6,k)),2)) for k in range(2,9)}
print('with122',OUT['with122'])
# harmonic baseline
t=np.arange(N); H2=np.column_stack([np.ones(N)]+[f(2*np.pi*h*t/365.25) for h in (1,2) for f in (np.sin,np.cos)])
H3=np.column_stack([H2]+[f(2*np.pi*3*t/365.25) for f in (np.sin,np.cos)])
fit2=(H2@np.linalg.lstsq(H2,Z,rcond=None)[0]).T; fit3=(H3@np.linalg.lstsq(H3,Z,rcond=None)[0]).T
corr=[np.corrcoef(main[c],fit2[c])[0,1] for c in range(C)]
OUT['harm_corr']=dict(min=round(min(corr),3),median=round(float(np.median(corr)),3))
mon=w.resample('MS').sum(); clim=mon.groupby(mon.index.month).mean().T.values
rawl={k:ward(Z.T,k) for k in range(2,9)}
tab=[]
for nmx,F in [('MSSA',main),('Harmonic2',fit2),('Harmonic3',fit3),('MSSA+122',F6),('Raw',Z.T),('Clim',clim)]:
    for k in range(2,9):
        l,s,ch=score(F,k)
        tab.append(dict(input=nmx,k=k,sil=round(s,3),ch=round(ch,1),sizes='–'.join(map(str,sorted(np.bincount(l)[1:],reverse=True))),ari=round(ari(base[k],l),2)))
OUT['table']=tab
for r in tab: print(r)
# transforms
for tr,fn in [('log1p',np.log1p),('sqrt',np.sqrt)]:
    s2,R2=mssa(stdz(fn(P)),365,8); ps2=[per(r.mean(1)) for r in R2]; sel=[j for j,p in enumerate(ps2) if 160<p<200 or 330<p<400]
    F=sum(R2[j] for j in sel).T
    OUT[tr]=dict(selected=[j+1 for j in sel],share=round(100*s2[sel].sum(),2),**{f'k{k}':dict(sil=round(score(F,k)[1],3),ari=round(ari(base[k],ward(F,k)),2)) for k in (3,4)})
    print(tr,OUT[tr])
json.dump(OUT,open('rev/a1.json','w'),indent=1)
pickle.dump(dict(base=base,fit2=fit2,clim=clim,main=main),open('rev/a1.pkl','wb'))
