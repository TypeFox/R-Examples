######################################################
# Section 5.6  The Example
######################################################

library(LearnBayes)

data(cancermortality)

 fit=laplace(betabinexch,c(-7,6),cancermortality)
 fit

npar=list(m=fit$mode,v=fit$var)
mycontour(lbinorm,c(-8,-4.5,3,16.5),npar,
  xlab="logit eta", ylab="log K")

 se=sqrt(diag(fit$var))
 fit$mode-1.645*se
 fit$mode+1.645*se

