import sys,pickle,json,collections;sys.path.insert(0,'/tmp/claude-0/-home-claude/6cffcb31-050c-5a22-af60-b2eb00fe5b65/scratchpad')
from rssa_like import *
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
_,_,names,w=pickle.load(open('chirps.pkl','rb')); names=list(names)
sens=pickle.load(open('sens.pkl','rb')); A=pickle.load(open('rev/a1.pkl','rb')); sh,R=pickle.load(open('rev/mssa365.pkl','rb'))
P=w.values;N,C=P.shape
def stdz(X): return (X-X.mean(0))/X.std(0,ddof=1)
def rowz(F): return (F-F.mean(1,keepdims=True))/F.std(1,ddof=1,keepdims=True)
def per(x): x=x-x.mean();f=np.abs(np.fft.rfft(x))[1:];return len(x)/(f.argmax()+1)
def ward(F,k,m='ward'): return fcluster(linkage(rowz(F),m),k,'maxclust')
def seasF(Zs):
    s2,R2=mssa(Zs,365,8); ps=[per(r.mean(1)) for r in R2]; g=[j for j,p in enumerate(ps) if 160<p<200 or 330<p<400]
    return sum(R2[j] for j in g).T
main=A['main'];fit2=A['fit2'];clim=A['clim'];Z=stdz(P);yrs=w.index.year.values
t=np.arange(N);H3=np.column_stack([np.ones(N)]+[f(2*np.pi*h*t/365.25) for h in (1,2,3) for f in (np.sin,np.cos)]);fit3=(H3@np.linalg.lstsq(H3,Z,rcond=None)[0]).T
k=4; base=ward(main,k); m={2:1,1:2,3:3,4:4}; base=np.array([m[b] for b in base])
V={'L180':sens[180][1],'L730':sens[730][1],'L1095':sens[1095][1],'avg':ward(main,k,'average'),'complete':ward(main,k,'complete'),
   'kmeans':KMeans(k,n_init=50,random_state=1).fit_predict(rowz(main))+1,'sub1':ward(seasF(stdz(P[yrs<=2009])),k),'sub2':ward(seasF(stdz(P[yrs>=2010])),k),
   '122':ward(sum(R[j] for j in [0,1,2,3,5,6]).T,k),'log':ward(seasF(stdz(np.log1p(P))),k),'sqrt':ward(seasF(stdz(np.sqrt(P))),k),'harm2':ward(fit2,k),'harm3':ward(fit3,k),'clim':ward(clim,k)}
alt=collections.defaultdict(list)
for nm,l in V.items():
    M=np.zeros((k,k))
    for a,b in zip(base,l): M[a-1,b-1]+=1
    r,c=linear_sum_assignment(-M); mp={c[i]+1:r[i]+1 for i in range(k)}
    al=[mp[x] for x in l]
    for n,b,a in zip(names,base,al):
        if a!=b: alt[n].append(a)
    if nm=='harm2': print('harm2 regions:',{g:[n for n,a in zip(names,al) if a==g] for g in range(1,5)})
for n in names:
    if alt[n]: print(n,'base',base[names.index(n)],'moves to',collections.Counter(alt[n]).most_common())
