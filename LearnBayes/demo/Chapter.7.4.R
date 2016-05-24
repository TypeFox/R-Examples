##############################################
# Section 7.4 Equal Mortality Rates?
##############################################

library(LearnBayes)
data(hearttransplants)
attach(hearttransplants)

sum(y)
sum(e)

 lambda=rgamma(1000,shape=277,rate=294681)
 ys94=rpois(1000,e[94]*lambda)

 hist(ys94,breaks=seq(0.5,max(ys94)+0.5))
 lines(c(y[94],y[94]),c(0,120),lwd=3)

S=readline(prompt="Type  <Return>   to continue : ")

lambda=rgamma(1000,shape=277,rate=294681)
prob.out=function(i)
{
   ysi=rpois(1000,e[i]*lambda)
   pleft=sum(ysi<=y[i])/1000
   pright=sum(ysi>=y[i])/1000
   min(pleft,pright)
 }
pout=sapply(1:94,prob.out)
windows()
plot(log(e),pout,ylab="Prob(extreme)")
