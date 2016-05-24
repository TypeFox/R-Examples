"UPMEpik2frompikw" <-function(pik,w)
{
n=sum(pik)
n=.as_int(n)
N=length(pik)
M=array(0,c(N,N))
for(k in 1:N) for(l in 1:N)
  if(pik[k]!=pik[l] & k!=l) M[k,l]= (pik[k]*w[l]-pik[l]*w[k])/(w[l]-w[k]) else M[k,l]=-1
for(i in 1:N) M[i,i]=pik[i]
for(k in 1:N)
  {
  tt=0
  comp=0
  for(l in 1:N)
        {if(M[k,l]!=-1) tt=tt+M[k,l]
         else comp=comp+1
        }
        cc=(n*pik[k]-tt)/comp
for(l in 1:N)  if(M[k,l]==-1) M[k,l]=cc
  }
M
}

