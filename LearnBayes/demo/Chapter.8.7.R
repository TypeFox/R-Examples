###################################################
# Section 8.7  Is a Baseball Hitter Really Streaky?
###################################################

library(LearnBayes)

 data(jeter2004)
 attach(jeter2004)
 data=cbind(H,AB)
 data1=regroup(data,5)

 log.marg=function(logK) 
     laplace(bfexch,0,list(data=data1,K=exp(logK)))$int
 
 log.K=seq(2,6)
 K=exp(log.K)
 log.BF=sapply(log.K,log.marg)
 BF=exp(log.BF)
 round(data.frame(log.K,K,log.BF,BF),2)