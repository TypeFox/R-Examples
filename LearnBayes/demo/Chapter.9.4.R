##############################################
# Section 9.4 Survival Modeling
##############################################

 library(LearnBayes)

 data(chemotherapy)
 attach(chemotherapy)
 library(survival)
 survreg(Surv(time,status)~factor(treat)+age,dist="weibull")

 start=c(-.5,9,.5,-.05)
 d=cbind(time,status,treat-1,age)
 fit=laplace(weibullregpost,start,d)
 fit

 proposal=list(var=fit$var,scale=1.5)
 bayesfit=rwmetrop(weibullregpost,proposal,fit$mode,10000,d)
 bayesfit$accept

 par(mfrow=c(2,2))
 sigma=exp(bayesfit$par[,1])
 mu=bayesfit$par[,2]
 beta1=bayesfit$par[,3]
 beta2=bayesfit$par[,4]
 hist(beta1,xlab="treatment",main="")
 hist(beta2,xlab="age",main="")
 hist(sigma,xlab="sigma",main="")
