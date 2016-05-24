gibbsmix=function(Niter=10^4,v=1){
#simulating the posterior of a mixture likelihood by Gibbs sampling

#data
da=sample(c(rnorm(10^2),2.5+rnorm(4*10^2)))

#log-likelihood function
like=function(mu){
  sum(log((.2*dnorm(da-mu[1])+.8*dnorm(da-mu[2]))))
  }

#log-likelihood surface
mu1=mu2=seq(-2,5,le=250)
lli=matrix(0,nco=250,nro=250)

for (i in 1:250)
for (j in 1:250)
  lli[i,j]=like(c(mu1[i],mu2[j]))

x=prop=runif(2,-2,5)
the=matrix(x,ncol=2)
curlike=hval=like(x)

for (i in 2:Niter){

# Gibbs samplin
  pp=1/(1+((.8*dnorm(da,mean=the[i-1,2]))/(.2*dnorm(da,mean=the[i-1,1]))))
  z=2-(runif(length(da))<pp)
  prop[1]=(v*sum(da[z==1]))/(sum(z==1)*v+1)+rnorm(1)*sqrt(v/(1+sum(z==1)*v))
  prop[2]=(v*sum(da[z==2]))/(sum(z==2)*v+1)+rnorm(1)*sqrt(v/(1+sum(z==2)*v))

  curlike=like(prop);hval=c(hval,curlike);the=rbind(the,prop)
  }

par(mar=c(4,4,1,1))
image(mu1,mu2,lli,xlab=expression(mu[1]),ylab=expression(mu[2]))
contour(mu1,mu2,-lli,nle=100,add=T)
points(the[,1],the[,2],cex=.6,pch=19)
lines(the[,1],the[,2],cex=.6,pch=19)
}
