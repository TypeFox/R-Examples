sunterpi2 <-
function(x,n)
{
xo<-sort(x,decreasing=TRUE)
if(xo[1]*n/sum(xo)>1) stop("There are some units with inclusion probability >1")
N<-length(xo)
t<-rev(cumsum(rev(xo)))
xbar<-t/(N:1)
kk0<-n*xo/t
k0<-which(kk0>=1,arr.ind = FALSE)[1]
kstar<-min(k0,N-n+1)
g<-numeric(kstar)
g[1]<-1/t[2]
for(i in 2:kstar)
{
g[i]<-prod(1-xo[1:(i-1)]/t[2:i])/t[i+1]
}

piij<-matrix(0,N,N)

for(k in 1:(N-1))
{
for(l in (k+1):N)
{
if(k<l & l<kstar)
{
piij[k,l]<-n*(n-1)/t[1]*g[k]*xo[k]*xo[l]
}
else if(k<kstar & kstar<=l)
{
piij[k,l]<-n*(n-1)/t[1]*g[k]*xo[k]*xbar[kstar]
}
else if(kstar<=k & k<l)
{
piij[k,l]<-n*(n-1)/t[1]*g[kstar-1]*(t[kstar]-xo[kstar-1])/(t[kstar]-xbar[kstar])*xbar[kstar]^2
}
piij[l,k]<-piij[k,l]
}
}
ord<-N+1-rank(x,ties.method="random")
piij<-piij[ord,ord]
diag(piij)<-apply(piij,1,sum)/(n-1)
piij
}

