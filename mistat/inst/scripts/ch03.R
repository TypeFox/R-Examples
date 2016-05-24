###################################################
### Chap03Start
###################################################



###################################################
### PlotDiscreteDistribution
###################################################
layout(matrix(1:2, nrow=1))
X <- 0:5
plot(X, dbinom(X, size=5, prob=0.5), 
     type="h", xlab="", ylab="")
plot(X, pbinom(X, size=5, prob=0.5), 
     type="s", xlab="", ylab="")
layout(1)
rm(X)


###################################################
### PlotContinuousDistribution
###################################################
X <- seq(0, 5, by=0.05)
plot(X, 1-exp(-X), 
     type="l", ylab=expression(F(x)), xlab=expression(x))
rm(X)


###################################################
### Binom
###################################################
X <- data.frame(i=0:30, 
                b=dbinom(x=0:30, size=30, prob=0.6),
                B=pbinom(q=0:30, size=30, prob = 0.6))
rm(X)


###################################################
### BinomInvCdf
###################################################
qbinom(p=0.5, size=30, prob=0.6)


###################################################
### PlotBinomialDistribution
###################################################
X <- 0:50
plot(X, dbinom(X, size=50, prob=0.25), 
     type="l", ylim=c(0, 0.14), 
     xlab=expression(i),
     ylab=expression(b(i, 50, p)))
lines(X, dbinom(X, size=50, prob=0.5), 
      type="l", lty="dashed")
lines(X, dbinom(X, size=50, prob=0.75), 
      type="l", lty="dotdash")
text(12, 0.135, expression(p==0.25))
text(25, 0.12, expression(p==0.5))
text(38, 0.135, expression(p==0.75))
rm(X)


###################################################
### Hyper
###################################################
X <- data.frame(j=0:10, 
                h=round(dhyper(x=0:10, m=15, n=60, k=10), 4),
                H=round(phyper(q=0:10, m=15, n=60, k=10), 4))
X
rm(X)


###################################################
### PlotHypergeometricDistribution
###################################################
X <- 0:100
plot(X, dhyper(X, m=500+350, n=350, k=100), 
     type="l", 
     xlab=expression(x), 
     ylab=expression("Hypergeometric p.d.f."))
rm(X)


###################################################
### PlotPoisoonDistribution
###################################################
X <- 0:50
plot(X, dpois(X, lambda=5), 
     type="l", ylim=c(0, 0.2), 
     xlab=expression(i),
     ylab=expression(p(i, lambda)))
lines(X, dpois(X, lambda=10), 
      type="l", lty="dashed")
lines(X, dpois(X, lambda=15), 
      type="l", lty="dotdash")
text(9, 0.19, expression(lambda==5))
text(14, 0.14, expression(lambda==10))
text(19, 0.12, expression(lambda==15))
rm(X)


###################################################
### PlotBinomialDistribution
###################################################
X <- 0:100
plot(X, dnbinom(X, size=5, prob=0.20), 
     type="l", 
     xlab=expression(i),
     ylab=expression(nb(i, 5, p)))
lines(X, dnbinom(X, size=5, prob=0.1), 
      type="l", lty="dashed")
text(35, 0.04, expression(p==0.20))
text(70, 0.015, expression(p==0.10))
rm(X)


###################################################
### PlotNormalDistribution
###################################################
X <- seq(0, 20, by=0.1)
plot(X, dnorm(X, mean=10, sd=1), 
     type="l", 
     xlab=expression(x),
     ylab=expression(f(x)))
lines(X, dnorm(X, mean=10, sd=2), 
      type="l", lty="dashed")
lines(X, dnorm(X, mean=10, sd=3), 
      type="l", lty="dotdash")
arrows(x1=c(11, 11, 11), 
       y1=c(
         dnorm(11, 10, 1),
         dnorm(11, 10, 2),
         dnorm(11, 10, 3)), 
       x0=c(16, 16, 16),
       y0=c(
         dnorm(11, 10, 1)+0.1,
         dnorm(11, 10, 2)+0.1,
         dnorm(11, 10, 3)+0.1))
text(c(18, 18, 18), 
     c(
       dnorm(11, 10, 1)+0.1,
       dnorm(11, 10, 2)+0.1,
       dnorm(11, 10, 3)+0.1), 
     c(expression(sigma==1), 
       expression(sigma==2), 
       expression(sigma==3)))
rm(X)


###################################################
### PlotStandardNormalDistribution
###################################################
X <- seq(-3, 3, by=0.1)
plot(X, pnorm(X, mean=1, sd=1), 
     type="l", 
     xlab=expression(x),
     ylab=expression("Standard Normal c.d.f."))
rm(X)


###################################################
### NormalProb
###################################################
pnorm(q=1.5, mean=0, sd=1)


###################################################
### NormQuant
###################################################
qnorm(p=0.95, mean=0, sd=1)


###################################################
### Example3_25
###################################################
pnorm(q=60.1, mean=60.02, sd=0.048, lower.tail=TRUE)


###################################################
### PlotLogNormIndexOfSkewness
###################################################
Sigma <- seq(from=sqrt(0.1), to=sqrt(3), length.out=200)

Skewness <- (exp(3*Sigma^2) -(3 * exp(Sigma^2)) + 2)/((exp(Sigma^2)-1)^(3/2)) 

plot(Sigma^2, Skewness, type="l", xlab=expression(sigma^2))

rm(Sigma, Skewness)


###################################################
### PlotExponentialDistribution
###################################################
X <- seq(0, 10, by=0.1)
plot(X, dexp(X, rate=1), 
     type="l", 
     xlab=expression(x),
     ylab=expression(f(x)))
lines(X, dexp(X, rate=1/2), 
      type="l", lty="dashed")
lines(X, dexp(X, rate=1/3), 
      type="l", lty="dotdash")
arrows(x1=c(0.5, 0.5, 0.5), 
       y1=c(
         dexp(0.5, rate=1),
         dexp(0.5, rate=1/2),
         dexp(0.5, rate=1/3)), 
       x0=c(4, 4, 4),
       y0=c(
         dexp(0.5, rate=1)+0.1,
         dexp(0.5, rate=1/2)+0.1,
         dexp(0.5, rate=1/3)+0.1))
text(c(6, 6, 6), 
     c(
       dexp(0.5, rate=1)+0.1,
       dexp(0.5, rate=1/2)+0.1,
       dexp(0.5, rate=1/3)+0.1), 
     c(expression(beta==1), 
       expression(beta==2), 
       expression(beta==3)))
rm(X)


###################################################
### GammaP
###################################################
pgamma(q=1, shape=1, scale=1)


###################################################
### GammaFunction
###################################################
gamma(5)


###################################################
### PlotGammaDistribution
###################################################
X <- seq(0, 6, by=0.1)
plot(X, dgamma(X, shape=0.5), 
     type="l", 
     xlab=expression(x),
     ylab=expression(f(x)))
lines(X, dgamma(X, shape=1), 
      type="l", lty="dashed")
lines(X, dgamma(X, shape=2), 
      type="l", lty="dotdash")
arrows(x1=c(1.5, 1.5, 1.5), 
       y1=c(
         dgamma(1.5, shape=0.5),
         dgamma(1.5, shape=1),
         dgamma(1.5, shape=2)), 
       x0=c(4, 4, 4),
       y0=c(
         dgamma(1.5, shape=0.5)+0.3,
         dgamma(1.5, shape=1)+0.3,
         dgamma(1.5, shape=2)+0.3))
text(c(5, 5, 5), 
     c(
       dgamma(1.5, shape=0.5)+0.3,
       dgamma(1.5, shape=1)+0.3,
       dgamma(1.5, shape=2)+0.3), 
     c(expression(beta==0.5), 
       expression(beta==1), 
       expression(beta==2)))
rm(X)


###################################################
### PlotWeibullDistribution
###################################################
X <- seq(0, 3, by=0.01)
plot(X, dweibull(X, shape=2), 
     type="l", 
     xlab=expression(x),
     ylab=expression(w(x,alpha,beta==1)))
lines(X, dweibull(X, shape=1.5), 
      type="l", lty="dashed")
arrows(x1=c(1, 1), 
       y1=c(
         dweibull(1, shape=2),
         dweibull(1, shape=1.5)), 
       x0=c(2, 2),
       y0=c(
         dweibull(1, shape=2)+0.1,
         dweibull(1, shape=1.5)+0.1))
text(c(2.5, 2.5), 
     c(
       dweibull(1, shape=2)+0.1,
       dweibull(1, shape=1.5)+0.1), 
     c(expression(alpha==2), 
       expression(alpha==1.5)))
rm(X)


###################################################
### PlotBetaDistribution
###################################################
X <- seq(0, 1, length.out=200)
plot(X, dbeta(X, shape1=2.5, shape2=5), 
     type="l", 
     xlab=expression(x),
     ylab=expression(w(x,v1,v2)))
lines(X, dbeta(X, shape1=2.5, shape2=2.5), 
      type="l", lty="dashed")
# arrows(x1=c(1, 1), 
#        y1=c(
#          dbeta(X, shape1=2.5, shape2=5),
#          dbeta(X, shape1=2.5, shape2=2.5)), 
#        x0=c(2, 2),
#        y0=c(
#          dbeta(X, shape1=2.5, shape2=5)+0.1,
#          dbeta(X, shape1=2.5, shape2=2.5)+0.1))
text(c(0.65, 0.75), 
     c(
       dbeta(0.4, shape1=2.5, shape2=5)+0.25,
       dbeta(0.6, shape1=2.5, shape2=2.5)+0.2), 
     c(expression("v1 = 2.5, v2 = 5"), 
       expression("v1 = 2.5, v2 = 2.5")))
rm(X)


###################################################
### PlotBivariateNormalDistribution
###################################################
library(mvtnorm)
X <- seq(-3, 3, by=0.25)
Y <- seq(-3, 3, by=0.25)
dmvnormOuter <- function(x,y){
  dmvnorm(cbind(x,y))
}
persp(X, Y, outer(X, Y, dmvnormOuter), 
      zlab=expression(density), 
      theta=50, phi=25, 
      axes=TRUE, ticktype="detailed")
rm(X, Y, dmvnormOuter)


###################################################
### PlotSeriesAndParallelSystem
###################################################
ParMar <- par("mar")

par(mar=c(0,0,0,0))

plot(1:9, 1:9, type="n", axes=FALSE, xlab="", ylab="")

rect(xleft=3, ybottom=7, xright=4, ytop=8)

rect(xleft=6, ybottom=7, xright=7, ytop=8)

rect(xleft=4.5, ybottom=1, xright=5.5, ytop=2)

rect(xleft=4.5, ybottom=3, xright=5.5, ytop=4)

arrows(1, 7.5, 3, 7.5)

arrows(4, 7.5, 6, 7.5)

arrows(7, 7.5, 9, 7.5)


arrows(1, 2.5, 3, 2.5)

arrows(3, 3.5, 4.5, 3.5)

arrows(3, 1.5, 4.5, 1.5)

arrows(7, 2.5, 9, 2.5)

lines(c(3, 3), c(1.5, 3.5))

lines(c(7, 7), c(1.5, 3.5))

lines(c(5.5, 7), c(1.5, 1.5))

lines(c(5.5, 7), c(3.5, 3.5))

text(5, 9, labels="Components in Series", cex=1.5)

text(5, 5, labels="Components in Parallel", cex=1.5)

text(3.5, 7.5, labels="C1", cex=1)

text(6.5, 7.5, labels="C2", cex=1)

text(5, 3.5, labels="C1", cex=1)

text(5, 1.5, labels="C2", cex=1)

invisible(par(ParMar))


###################################################
### PlotStudentTDistribution
###################################################
X <- seq(-5, 5, by=0.1)
plot(X, dt(X, df=5), 
     type="l", 
     xlab=expression(x),
     ylab=expression(t(x,v)),
     ylim=c(0, 0.5))
lines(X, dt(X, df=50), 
      type="l", lty="dashed")

arrows(x1=c(0, 2.5), 
       y1=c(
         dt(0, df=50),
         dt(2.5, df=5)), 
       x0=c(3, 3),
       y0=c(
         dt(0, df=50)+0.05,
         dt(2.5, df=5)+0.05))
text(c(4, 4), 
     c(
       dt(0, df=50)+0.06,
       dt(2.5, df=5)+0.06), 
     c(expression(v==50), 
       expression(v==5)))
rm(X)


###################################################
### PlotFDistribution
###################################################
X <- seq(0, 6, by=0.05)
plot(X, df(X, df1=10, df2=10), 
     type="l", 
     xlab=expression(x),
     ylab=expression("F(v1,v2)"))
rm(X)


###################################################
### Chap03End
###################################################
detach(package:mvtnorm)
