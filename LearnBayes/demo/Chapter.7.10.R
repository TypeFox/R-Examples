#################################################
# Section 7.10 Posterior Predictive Model Checking
#################################################

library(LearnBayes)
 data(hearttransplants)
 attach(hearttransplants)

 datapar = list(data = hearttransplants, z0 = 0.53)

 start = c(4, -7)
 fitgibbs = gibbs(poissgamexch, start, 1000, c(1,.15), datapar)

 lam94=rgamma(1000,y[94]+alpha,e[94]+alpha/mu)

 ys94=rpois(1000,e[94]*lam94)

 hist(ys94,breaks=seq(-0.5,max(ys94)+0.5))
 lines(y[94]*c(1,1),c(0,100),lwd=3)

S=readline(prompt="Type  <Return>   to continue : ")

prob.out=function(i)
{
   lami=rgamma(1000,y[i]+alpha,e[i]+alpha/mu)
   ysi=rpois(1000,e[i]*lami)
   pleft=sum(ysi<=y[i])/1000
   pright=sum(ysi>=y[i])/1000
   min(pleft,pright)
 }
pout.exchange=sapply(1:94,prob.out)

 windows()
 plot(pout,pout.exchange,xlab="P(extreme), equal means",
 ylab="P(extreme), exchangeable")
 abline(0,1)


