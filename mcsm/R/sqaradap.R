sqaradap=function(T=10^4,TT=10^4,scale=.5,factor=1){
#An adaptive example on the noisy squared AR
library(coda)

rho=.85
tau=.2

#ensure the posterior is bimodal and not centered near sqrt(yc)
xc=1
while (xc>0){

  xm=rnorm(1,mean=-2)
  xc=rho*xm+rnorm(1);xp=rho*xc+rnorm(1)
  yc=xc^2+tau*rnorm(1)
  }

#targets, log and not log:
ef=function(x){

    -.5*((xm*rho-x)^2+(x*rho-xp)^2+(yc-x^2)^2/tau^2 )
    }

eef=function(x){exp(ef(x))}

#Initial sample
xmc=sqrt(abs(yc))
for (t in 2:T){

     prop=xmc[t-1]+scale*rnorm(1)
     if (log(runif(1))>ef(prop)-ef(xmc[t-1]))
        prop=xmc[t-1]

     xmc=c(xmc,prop)
     }

#Use previous samples to build NP proposal
bw=factor*bw.nrd0(xmc)
curdens=log(density(xmc,from=xmc[T],to=xmc[T],n=1,bw=bw)$y)

for (t in (T+1):TT){

  bw=500*bw.nrd0(xmc)
  prop=rnorm(1,mean=sample(xmc,1),sd=bw)
  prodens=log(density(xmc,from=prop,to=prop,n=1,bw=bw)$y)
  if ((is.na(prop))||(prodens==-Inf)||(log(runif(1))>ef(prop)-ef(xmc[t-1])+curdens-prodens)){

     prop=xmc[t-1];prodens=curdens}

  xmc=c(xmc,prop)
  curdens=prodens
} 

mymc=xmc[(TT/10):TT]
par(mfrow=c(1,2),mar=c(4,2,1,1))
plot(xmc[seq(1,length(xmc),le=5000)],cex=.2,type="l",xlab="Iterations x 20")
hist(mymc,nclass=100,col="grey89",pro=T,xlab=expression(x[t]),ylab="",main="")
ordin=apply(as.matrix(seq(min(mymc),max(mymc),le=200)),1,eef)
lines(seq(min(mymc),max(mymc),le=200),ordin*max(density(mymc)$y)/max(ordin),lwd=2,col="gold4")
}
