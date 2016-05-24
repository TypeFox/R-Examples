##############################################
# Section 5.10 Sampling Importance Resampling
##############################################

library(LearnBayes)
data(cancermortality)
fit=laplace(betabinexch,c(-7,6),cancermortality)

tpar=list(m=fit$mode,var=2*fit$var,df=4)

theta.s=sir(betabinexch,tpar,10000,cancermortality)

 S=bayes.influence(theta.s,cancermortality)

 plot(c(0,0,0),S$summary,type="b",lwd=3,xlim=c(-1,21),
  ylim=c(5,11), xlab="Observation removed",ylab="log K")
 for (i in 1:20)
  lines(c(i,i,i),S$summary.obs[i,],type="b")
