#################################################
# Section 8.6 Models for Soccer Goals
#################################################

library(LearnBayes)

 data(soccergoals)
 attach(soccergoals)
 datapar=list(data=goals,par=c(4.57,1.43))
 fit1=laplace(logpoissgamma,.5,datapar)
 datapar=list(data=goals,par=c(1,.5))
 fit2=laplace(logpoissnormal,.5,datapar)
 datapar=list(data=goals,par=c(2,.5))
 fit3=laplace(logpoissnormal,.5,datapar)
 datapar=list(data=goals,par=c(1,2))
 fit4=laplace(logpoissnormal,.5,datapar)

 postmode=c(fit1$mode,fit2$mode,fit3$mode,fit4$mode)
 postsd=sqrt(c(fit1$var,fit2$var,fit3$var,fit4$var))
 logmarg=c(fit1$int,fit2$int,fit3$int,fit4$int)
 cbind(postmode,postsd,logmarg)

