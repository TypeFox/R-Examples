estN.bootstrap <-
function(z, rnhat, nboot = 200){
rnhat=round(rnhat)
m=length(z)
zcons=rep(0,m+1)
zcons[m+1]=rnhat-sum(z)
for(i in 1:m){zcons[i]=zcons[m+1]+sum(z[1:i])}
zcons=zcons/rnhat

zstar=matrix(0,m,nboot)
for(j in 1:nboot){
a=runif(rnhat,0,1)
for(i in 1:m){
b=sort(a)
if(i==1) between=b[(b>=zcons[m+1]) & (b<zcons[1])] else between=b[(b>=zcons[i-1]) & (b<zcons[i])]
zstar[i,j]=length(between)
}
}
zstar
}
