import sys,pickle,json;sys.path.insert(0,'/tmp/claude-0/-home-claude/6cffcb31-050c-5a22-af60-b2eb00fe5b65/scratchpad')
from rssa_like import *
from sklearn.metrics import adjusted_rand_score as ari
from sklearn.cluster import KMeans
from scipy.optimize import linear_sum_assignment
_,_,names,w=pickle.load(open('chirps.pkl','rb')); names=list(names)
sens=pickle.load(open('sens.pkl','rb')); sh,R=pickle.load(open('rev/mssa365.pkl','rb')); A=pickle.load(open('rev/a1.pkl','rb'))
P=w.values;N,C=P.shape
def stdz(X): return (X-X.mean(0))/X.std(0,ddof=1)
def rowz(F): return (F-F.mean(1,keepdims=True))/F.std(1,ddof=1,keepdims=True)
def per(x): x=x-x.mean();f=np.abs(np.fft.rfft(x))[1:];return len(x)/(f.argmax()+1)
def ward(F,k,m='ward'): return fcluster(linkage(rowz(F),m),k,'maxclust')
def seasF(Zs,L=365):
    s2,R2=mssa(Zs,L,8); ps=[per(r.mean(1)) for r in R2]; g=[j for j,p in enumerate(ps) if 160<p<200 or 330<p<400]
    return sum(R2[j] for j in g).T
main=A['main']; fit2=A['fit2']; clim=A['clim']
t=np.arange(N); H3=np.column_stack([np.ones(N)]+[f(2*np.pi*h*t/365.25) for h in (1,2,3) for f in (np.sin,np.cos)])
Z=stdz(P); fit3=(H3@np.linalg.lstsq(H3,Z,rcond=None)[0]).T
ps=[per(r.mean(1)) for r in R]; F122=sum(R[j] for j in [0,1,2,3,5,6]).T
Flog=seasF(stdz(np.log1p(P))); Fsq=seasF(stdz(np.sqrt(P)))
yrs=w.index.year.values
F1=seasF(stdz(P[yrs<=2009])); F2=seasF(stdz(P[yrs>=2010]))
# fold-mean of MSSA seasonal signal by day of year
doy=np.minimum(w.index.dayofyear.values,365)-1
fold=np.stack([np.bincount(doy,weights=main[c],minlength=365)/np.bincount(doy,minlength=365) for c in range(C)])
corr=[np.corrcoef(fold[c],(H3[:365,:5]@np.linalg.lstsq(H3[:,:5],Z[:,c],rcond=None)[0]))[0,1] for c in range(C)]
print('fold-mean vs 2-harmonic corr min/median',round(min(corr),4),round(float(np.median(corr)),4))
res={}
for k in (3,4):
    base=ward(main,k)
    V={'L = 180 days':sens[180][0 if k==3 else 1],'L = 730 days':sens[730][0 if k==3 else 1],'L = 1,095 days':sens[1095][0 if k==3 else 1],
       'Average linkage':ward(main,k,'average'),'Complete linkage':ward(main,k,'complete'),
       'k-means':KMeans(k,n_init=50,random_state=1).fit_predict(rowz(main))+1,
       'Sub-period 1995–2009':ward(F1,k),'Sub-period 2010–2024':ward(F2,k),
       'Adding the 122-day pair':ward(F122,k),'log(1 + x) transform':ward(Flog,k),'Square-root transform':ward(Fsq,k),
       'Harmonic regression (2 harmonics)':ward(fit2,k),'Harmonic regression (3 harmonics)':ward(fit3,k),
       'Monthly climatology':ward(clim,k),'MSSA fold-mean (one-year cycle)':ward(fold,k)}
    keep=np.zeros(C); rows=[]; co=np.zeros((C,C))
    for nm,l in V.items():
        M=np.zeros((k,k))
        for a,b in zip(base,l): M[a-1,b-1]+=1
        r,c=linear_sum_assignment(-M); mp={c[i]+1:r[i]+1 for i in range(k)}
        al=np.array([mp[x] for x in l]); keep+=(al==base); rows.append(dict(cfg=nm,ari=round(ari(base,l),2),moved=int((al!=base).sum())))
        co+=(l[:,None]==l[None,:])
    co/=len(V); stab=keep/len(V)
    from scipy.spatial.distance import squareform
    cons=fcluster(linkage(squareform(1-co,checks=False),'average'),k,'maxclust')
    res[k]=dict(rows=rows,mean_ari=round(np.mean([r['ari'] for r in rows]),2),median_ari=round(float(np.median([r['ari'] for r in rows])),2),
               stab={n:round(s,2) for n,s in zip(names,stab)},base={n:int(b) for n,b in zip(names,base)},
               ari_consensus=round(ari(base,cons),2),n_stable75=int((stab>=0.75).sum()),n_full=int((stab==1).sum()))
    print('k',k,'meanARI',res[k]['mean_ari'],'median',res[k]['median_ari'],'stable>=.75',res[k]['n_stable75'],'full',res[k]['n_full'],'consensusARI',res[k]['ari_consensus'])
    for r in rows: print('  ',r)
    for z in sorted(set(base)): print('  zone',z,[(n,round(s,2)) for n,b,s in zip(names,base,stab) if b==z])
# Hennig bootstrap Jaccard (resample units)
rng=np.random.default_rng(2026)
F=rowz(main)
for k in (3,4):
    base=ward(main,k); J=[[] for _ in range(k)]
    for b in range(500):
        idx=rng.choice(C,C,replace=True); u=np.unique(idx)
        lb=fcluster(linkage(F[idx],'ward'),k,'maxclust')
        lab_u={}
        for i,l in zip(idx,lb): lab_u[i]=l
        for g in range(1,k+1):
            A_=set(i for i in u if base[i]==g)
            if not A_: continue
            best=max(len(A_&set(i for i in u if lab_u[i]==h))/len(A_|set(i for i in u if lab_u[i]==h)) for h in range(1,k+1))
            J[g-1].append(best)
    res[k]['jaccard']={g+1:round(float(np.mean(J[g])),2) for g in range(k)}
    print('k',k,'bootstrap Jaccard',res[k]['jaccard'],'sizes',np.bincount(base)[1:])
json.dump(res,open('rev/a2.json','w'),indent=1,ensure_ascii=False)
