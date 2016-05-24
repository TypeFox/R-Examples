################################
# Section 2.4 Using a Beta Prior
#############################

library(LearnBayes)

quantile2=list(p=.9,x=.5)
quantile1=list(p=.5,x=.3)
ab=beta.select(quantile1,quantile2)

 a = ab[1]
 b = ab[2]
 s = 11
 f = 16
 curve(dbeta(x,a+s,b+f), from=0, to=1, 
       xlab="p",ylab="Density",lty=1,lwd=4)
 curve(dbeta(x,s+1,f+1),add=TRUE,lty=2,lwd=4)
 curve(dbeta(x,a,b),add=TRUE,lty=3,lwd=4)
 legend(.7,4,c("Prior","Likelihood","Posterior"),
     lty=c(3,2,1),lwd=c(3,3,3))

 1 - pbeta(0.5, a + s, b + f)

 qbeta(c(0.05, 0.95), a + s, b + f)

 ps = rbeta(1000, a + s, b + f)

 windows()
 hist(ps,xlab="p")

 sum(ps >= 0.5)/1000

 quantile(ps, c(0.05, 0.95))
