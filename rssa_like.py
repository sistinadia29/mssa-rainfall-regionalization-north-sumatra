import numpy as np, pandas as pd, sys
from numpy.lib.stride_tricks import sliding_window_view as swv
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import silhouette_score, calinski_harabasz_score
def mssa(Z,L,nk=8):
    N,C=Z.shape;K=N-L+1
    S=np.zeros((L,L))
    Xs=[swv(Z[:,c],L).T for c in range(C)]
    for X in Xs: S+=X@X.T
    ev,U=np.linalg.eigh(S); o=ev.argsort()[::-1]; ev,U=ev[o],U[:,o]
    def hank(M):
        out=np.zeros(N);cnt=np.zeros(N)
        for i in range(L): out[i:i+K]+=M[i];cnt[i:i+K]+=1
        return out/cnt
    R=[]
    for j in range(nk):
        u=U[:,[j]]; R.append(np.stack([hank(u@(u.T@X)) for X in Xs],1))
    return ev/ev.sum(),R
def evalk(F,lab):
    F=(F-F.mean(1,keepdims=True))/F.std(1,ddof=1,keepdims=True);Lk=linkage(F,'ward');res={}
    for k in range(2,9):
        l=fcluster(Lk,k,'maxclust');res[k]=l
        print(lab,k,round(silhouette_score(F,l),3),round(calinski_harabasz_score(F,l),1),'-'.join(map(str,sorted(np.bincount(l)[1:],reverse=True))))
    return res
def period(x):
    x=x-x.mean();f=np.abs(np.fft.rfft(x))[1:];return round(len(x)/(f.argmax()+1),1)
