###############################
# Section 10.2 Robust Modeling
###############################

library(LearnBayes)

 data(darwin)
 attach(darwin)
 fit=robustt(difference,4,10000)

 plot(density(fit$mu),xlab="mu")

 mean.lambda=apply(fit$lam,2,mean)
 lam5=apply(fit$lam,2,quantile,.05)
 lam95=apply(fit$lam,2,quantile,.95)

S=readline(prompt="Type  <Return>   to continue : ")

 windows()
 plot(difference,mean.lambda,lwd=2,ylim=c(0,3),ylab="Lambda")
 for (i in 1:length(difference))
   lines(c(1,1)*difference[i],c(lam5[i],lam95[i]))
 points(difference,0*difference-.05,pch=19,cex=2)
