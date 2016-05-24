#############################################################
# Section 6.10 Analysis of the Stanford Heart Transplant Data
#############################################################

library(LearnBayes)

data(stanfordheart)

 start=c(0,3,-1)
 laplacefit=laplace(transplantpost,start,stanfordheart)
 laplacefit

 proposal=list(var=laplacefit$var,scale=2)
 s=rwmetrop(transplantpost,proposal,start,10000,stanfordheart)
 s$accept

 par(mfrow=c(2,2))
 tau=exp(s$par[,1])
 plot(density(tau),main="TAU")
 lambda=exp(s$par[,2])
 plot(density(lambda),main="LAMBDA")
 p=exp(s$par[,3])
 plot(density(p),main="P")

 apply(exp(s$par),2,quantile,c(.05,.5,.95))

S=readline(prompt="Type  <Return>   to continue : ")

 par(mfrow=c(1,1))
  t=seq(1,240)
 p5=0*t; p50=0*t; p95=0*t
 for (j in 1:240)
 { S=(lambda/(lambda+t[j]))^p
   q=quantile(S,c(.05,.5,.95))
   p5[j]=q[1]; p50[j]=q[2]; p95[j]=q[3]}
 windows()
 plot(t,p50,type="l",ylim=c(0,1),ylab="Prob(Survival)",
   xlab="time")
 lines(t,p5,lty=2)
 lines(t,p95,lty=2)
