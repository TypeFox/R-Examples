########################################################
# Section 7.5 Modeling a Prior Belief of Exchangeability
########################################################

library(LearnBayes)

pgexchprior=function(lambda,pars)
{
alpha=pars[1]; a=pars[2]; b=pars[3]
(alpha-1)*log(prod(lambda))-(2*alpha+a)*log(alpha*sum(lambda)+b)
}

alpha=c(5,20,80,400); par(mfrow=c(2,2))
for (j in 1:4)
    mycontour(pgexchprior,c(.001,5,.001,5),c(alpha[j],10,10),
      main=paste("ALPHA = ",alpha[j]),xlab="LAMBDA 1",ylab="LAMBDA 2")


