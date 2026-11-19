import sys,pickle,json;sys.path.insert(0,'/tmp/claude-0/-home-claude/6cffcb31-050c-5a22-af60-b2eb00fe5b65/scratchpad')
from rssa_like import *
from sklearn.metrics import adjusted_rand_score as ari
_,_,names,w=pickle.load(open('chirps.pkl','rb')); names=list(names)
sens=pickle.load(open('sens.pkl','rb')); A=pickle.load(open('rev/a1.pkl','rb')); sh,R=pickle.load(open('rev/mssa365.pkl','rb'))
def rowz(F): return (F-F.mean(1,keepdims=True))/F.std(1,ddof=1,keepdims=True)
def ward(F,k): return fcluster(linkage(rowz(F),'ward'),k,'maxclust')
main=A['main'];fit2=A['fit2'];clim=A['clim']
V={'MSSA365':ward(main,4),'Clim':ward(clim,4),'Harm2':ward(fit2,4),'L180':sens[180][1],'L730':sens[730][1],'L1095':sens[1095][1]}
ks=list(V)
print('      '+' '.join(f'{k:>8}' for k in ks))
for a in ks: print(f'{a:>8}'+' '.join(f'{ari(V[a],V[b]):8.2f}' for b in ks))
# harmonic amplitude ratio & phase per unit, and temporal variability of MSSA seasonal amplitude
P=w.values;N=len(P);t=np.arange(N);Z=(P-P.mean(0))/P.std(0,ddof=1)
H=np.column_stack([np.ones(N)]+[f(2*np.pi*h*t/365.25) for h in (1,2) for f in (np.sin,np.cos)])
B=np.linalg.lstsq(H,Z,rcond=None)[0]
a1=np.hypot(B[1],B[2]);a2=np.hypot(B[3],B[4])
yrs=w.index.year.values; amp=np.array([[main[c][yrs==y].std() for y in range(1995,2025)] for c in range(len(names))])
cv=amp.std(1)/amp.mean(1)
for z in (1,2,3,4): pass
lab=V['MSSA365']; 
for g in sorted(set(lab)):
    i=lab==g; print(g,i.sum(),'a2/a1',round(float(np.mean(a2[i]/a1[i])),2),'CV of yearly seasonal amplitude',round(float(cv[i].mean()),2))
print('where Labuhanbatu goes at k=3 in variants: see a2')
