gibbs=function(niter,datha,mix){

n=length(datha)
meanp=mean(datha)
varp=var(datha)
z=rep(0,n)
ssiz=rep(0,mix$k)
nxj=ssiz
nxj2=ssiz
ssum=ssiz
mug=matrix(0,nrow=niter,ncol=mix$k)
sigg=mug
prog=mug

lolik=rep(0,niter)
lik=matrix(0,n,mix$k)

perms=matrix(unlist(permn(mix$k)),ncol=mix$k,byrow=T)
nperms=dim(perms)[1]
chibdeno=0

for (i in 1:niter){
for (t in 1:n){

  prob=mix$p*dnorm(datha[t],mean=mix$mu,sd=sqrt(mix$sig))
  if (sum(prob)==0) prob=rep(1,mix$k)
  z[t]=sample(x=1:mix$k,size=1,prob=prob)
  }
  for (j in 1:mix$k){

    ssiz[j]=sum(z==j)
    nxj[j]=sum(as.numeric(z==j)*datha)
    nxj2[j]=sum(as.numeric(z==j)*datha^2)
    }
  mug[i,]=rnorm(mix$k,mean=(meanp+nxj)/(1+ssiz),sd=1/sqrt((1+ssiz)/mix$sig))
  mix$mu=mug[i,]
  for (j in 1:mix$k) ssum[j]=sum(as.numeric(z==j)*(datha-mix$mu[j])^2)
  sigg[i,]=1/rgamma(mix$k,shape=10+0.5*ssiz,rate=varp+0.5*(mix$mu-meanp)^2+0.5*ssum)
  mix$sig=sigg[i,]
  prog[i,]=rdirichlet(1,par=ssiz+0.5)
  mix$p=prog[i,]
  for (j in 1:mix$k) lik[,j]=mix$p[j]*dnorm(x=datha,mean=mix$mu[j],sd=sqrt(mix$sig[j])) 
  lolik[i]=sum(log(apply(lik,1,sum)))

  for (j in 1:nperms)
    chibdeno=chibdeno+exp(sum(dnorm(mix$mu[perms[j,]],
  mean=(meanp+nxj)/(1+ssiz),sd=1/sqrt((1+ssiz)/mix$sig[perms[j,]]),log=TRUE))+
  sum(dgamma(1/mix$sig[perms[j,]],10+ssiz/2+0.5,varp^2+meanp^2/2+nxj2/2-
  (nxj+meanp)^2/(2*(ssiz+1)),log=TRUE)-2*log(mix$sig[perms[j,]]))+
  sum((ssiz+0.5-1)*log(mix$p[perms[j,]]))+lgamma(sum(ssiz+0.5))-sum(lgamma(ssiz+0.5)))
  }
chibdeno=chibdeno/(niter*nperms)
list(k=mix$k,mu=mug,sig=sigg,p=prog,lolik=lolik,deno=chibdeno)
}
