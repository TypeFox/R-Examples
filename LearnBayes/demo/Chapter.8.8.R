###################################################################
# Section 8.8 A Test of Independence in a Two-Way Contingency Table
###################################################################

library(LearnBayes)

 data=matrix(c(11,9,68,23,3,5),c(2,3))
 data

 chisq.test(data)

 a=matrix(rep(1,6),c(2,3))
 a

 ctable(data,a)

  log.K=seq(2,7)
  compute.log.BF=function(log.K)
     log(bfindep(data,exp(log.K),100000)$bf)
  log.BF=sapply(log.K,compute.log.BF)
  BF=exp(log.BF)

round(data.frame(log.K,log.BF,BF),2)