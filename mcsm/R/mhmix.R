mhmix=function(Niter=10^4,lange=FALSE,scale=1){
#simulating the posterior of a mixture likelihood

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

iter=1
x=runif(2,-2,5)
the=matrix(x,ncol=2)
curlike=hval=like(x)

if (!lange){# random walk

while (iter<Niter){

  prop=the[iter,]+rnorm(2)*scale
  if ((max(-prop)>2)||(max(prop)>5)||(log(runif(1))>like(prop)-curlike))
      prop=the[iter,]
      
  curlike=like(prop);hval=c(hval,curlike);the=rbind(the,prop)
  iter=iter+1
  }

par(mar=c(4,4,1,1))
image(mu1,mu2,lli,xlab=expression(mu[1]),ylab=expression(mu[2]))
contour(mu1,mu2,-lli,nle=100,add=T)
points(unique(the[,1]),unique(the[,2]),cex=.6,pch=19)

}else{#Langevin version

#log-likelihood function
gradlike=function(mu){
  gr=sum(.2*(da-mu[1])*dnorm(da-mu[1])/((.2*dnorm(da-mu[1])+.8*dnorm(da-mu[2]))))
  gr=c(gr,sum(.8*(da-mu[2])*dnorm(da-mu[2])/((.2*dnorm(da-mu[1])+.8*dnorm(da-mu[2])))))

  return(gr)
  }

curmean=x+.5*scale^2*gradlike(x)

# random walk with drift
while (iter<Niter){

  prop=curmean+rnorm(2)*scale
  meanprop=prop+.5*scale^2*gradlike(prop)
  if ((max(-prop)>2)||(max(prop)>5)||(log(runif(1))>like(prop)-curlike-
	sum(dnorm(prop,curmean,lo=T))+sum(dnorm(the[iter,],meanprop,lo=T)))){
      prop=the[iter,];meanprop=curmean}

  curlike=like(prop);curmean=meanprop;
  hval=c(hval,curlike);the=rbind(the,prop)
  iter=iter+1
  }

par(mar=c(4,4,1,1))
image(mu1,mu2,lli,xlab=expression(mu[1]),ylab=expression(mu[2]))
contour(mu1,mu2,lli,nle=100,add=T)
points(unique(the[,1]),unique(the[,2]),cex=.6,pch=19)
lines(unique(the[,1]),unique(the[,2]),cex=.6,pch=19)
}

}
