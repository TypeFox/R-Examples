#########################################################
# Section 5.8 Rejection Sampling
#########################################################

library(LearnBayes)

data(cancermortality)
fit=laplace(betabinexch,c(-7,6),cancermortality)

betabinT=function(theta,datapar)
{
data=datapar$data
tpar=datapar$par
d=betabinexch(theta,data)-dmt(theta,mean=c(tpar$m),
  S=tpar$var,df=tpar$df,log=TRUE)
return(d)
}

tpar=list(m=fit$mode,var=2*fit$var,df=4)
datapar=list(data=cancermortality,par=tpar)

 start=c(-6.9,12.4)
 fit1=laplace(betabinT,start,datapar)
 fit1$mode

 betabinT(fit1$mode,datapar)

 theta=rejectsampling(betabinexch,tpar,-569.2813,10000,cancermortality)
 dim(theta)

 mycontour(betabinexch,c(-8,-4.5,3,16.5),cancermortality,
   xlab="logit eta",ylab="log K")
 points(theta[,1],theta[,2])
