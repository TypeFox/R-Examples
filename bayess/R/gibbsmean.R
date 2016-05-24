gibbsmean=function(p,datha,niter=10^4)
{
n=length(datha)
z=rep(0,n)
ssiz=nxj=rep(0,2)
mug=matrix(0,nrow=niter+1,ncol=2)
mug[1,]=rep(mean(datha),2)
for (i in 2:(niter+1))
{
for (t in 1:n)
{
prob=c(p,1-p)*dnorm(datha[t],mean=mug[i-1,])
z[t]=sample(c(1,2),size=1,prob=prob)   
}
for (j in 1:2)
{
ssiz[j]=1+sum(z==j)
nxj[j]=sum(as.numeric(z==j)*datha)
}
mug[i,]=rnorm(2,mean=nxj/ssiz,sd=sqrt(1/ssiz))
}
mug
}
