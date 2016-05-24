"UPmaxentropypi2" <-function(pik)
{
n=sum(pik)
n=.as_int(n)
N=length(pik)
M=array(0,c(N,N))
if(n>=2)
{
pik2=pik[pik>0 & pik<1]
pikt=UPMEpiktildefrompik(pik2)
w=pikt/(1-pikt)
M[pik>0 & pik<1,pik>0 & pik<1]=UPMEpik2frompikw(pik2,w)
M[,pik==1]=pik
for(k in 1:N) if(pik[k]==1)  M[k,]=pik
}
if(n==1) for(k in 1:N) M[k,k]=pik[k]
M
}

