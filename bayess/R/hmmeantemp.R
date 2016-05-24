hmmeantemp=function(dat,niter,var=1,alpha=1)
{
lpost=function(x,mu,p=0.7,delta=0,lambda=1)
{
sum(log(p * dnorm(x, mu[1]) + (1 - p) * dnorm(x, mu[2]))) +
sum(log(dnorm(mu, delta, 1/sqrt(lambda))))
}

mu=matrix(0,niter,2)
mu[1,]=c(1,3)
for (i in 2:niter)
{
muprop=rnorm(2,mu[i-1,],sqrt(var))
bound=lpost(dat,muprop)-lpost(dat,mu[i-1,])
if (runif(1)<=exp(alpha*bound)) mu[i,]=muprop else
mu[i,]=mu[i-1,]
}
mu
}
