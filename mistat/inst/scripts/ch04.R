###################################################
### Chap04Start
###################################################
library(mistat)
library(boot)


###################################################
### PlotHistogram100SampleMeans
###################################################
set.seed(123)

X <- rnorm(100, 50, 25)

set.seed(123)                          

B <- boot(X, 
          statistic=function(x, i, n) mean(x[i[1:n]]), 
          R=100, n=10)$t

hist(B, 
     xlab=expression(bar(x)),
     main="")

rm(B)                                  # Keep X for next


###################################################
### PlotHistogram100SampleVariance
###################################################
set.seed(123)                          

B <- boot(X, 
          statistic=function(x, i, n) var(x[i[1:n]]), 
          R=100, n=10)$t

hist(B, 
     xlab=expression(italic(S)^2),
     main="")

rm(X, B)


###################################################
### PlotOCcurveBinomial
###################################################
library(AcceptanceSampling)

X <- OC2c(n=20, c= 4, type="binomial")

plot(X@pd, X@paccept, type = "l", 
     xlab = expression(p), 
     ylab = expression("OC(p)"), ylim = c(0,1))

rm(X)


###################################################
### PlotCriticalRegionsOneSidedZtest
###################################################
X <- seq(-3, 3, length.out=400)

plot(X, dnorm(X), type="l", axes=FALSE, xlab="", ylab="", 
     xlim=c(-4, 5), ylim=c(-0.2, 0.45))
lines(c(0,0), c(0, 1), lty="dashed")
lines(rep(qnorm(p=0.90), 2), c(0, 1), lty="dashed")

arrows(-4, 0, 4, 0, length = 0.10)

arrows(-4, -0.1, 4, -0.1, length = 0.10)
lines(c(0, 0), c(-0.115, -0.085))
lines(rep(qnorm(p=0.90), 2), c(-0.115, -0.085))

polygon(x=c(X[X>=qnorm(p=0.90)], 3, qnorm(p=0.90)), 
        y=c(dnorm(X[X>=qnorm(p=0.90)]), 0, 0), 
        col="grey80", border=NA)

text(4.5, 0, labels=expression(paste(" ", italic(z), " scale")))
text(4.5, -0.1, labels=expression(paste(" ", italic(X), " scale")))

text(0, -0.15, labels=expression(mu[0]))
text(qnorm(p=0.90), -0.15, labels=expression(mu[0]+z[1-alpha]*frac(sigma, sqrt(n))))

text(0, -0.05, labels=expression(0))
text(qnorm(p=0.90), -0.05, labels=expression(z[1-alpha]))
text(-1.5, -0.05, labels="Accept")
text(2.5, -0.05, labels="Reject")

text(-2, 0.4, labels=expression(paste("acceptance\nprobability ", (1-alpha))))
arrows(-2, 0.35, -1, 0.1, length = 0.10)
text(3, 0.2, labels=expression(paste("rejection\nprobability ", alpha)))
arrows(3, 0.15, 1.8, 0.025, length = 0.10)

rm(X)


###################################################
### InvCdfBinom
###################################################
## qbinom(p=0.95, size=20, prob=0.2)


###################################################
### PlotSimulatedConfidenceIntervalsMean
###################################################
set.seed(123)

X <- rnorm(500, 10, 1)

library(boot)

confintFunc <- function(x, i, n, level=0.95){
  as.numeric(confint(lm(x[i[1:n]]~1), level=level))
}

set.seed(123)                          

B <- boot(X, statistic=confintFunc, R=50, n=10)$t

plot(c(1,50), c(min(B), max(B)), type = "n", 
     xlab = "", 
     ylab = "")

B <- cbind(B, 1:50)

B <- cbind(B, 1:50)

tmpFun <- function(k) lines(k[3:4], k[1:2])

invisible(apply(B, 1, tmpFun))

abline(h=10, lty="dashed")

rm(B, X, confintFunc, tmpFun)


###################################################
### PlotSimulated50qqnorm
###################################################
set.seed(123)

X <- rnorm(50, 10, 1)

qqnorm(y=X)

abline(a=10, b=1)

rm(X)


###################################################
### QqPlot
###################################################
library(car)

set.seed(123)

X <- rlnorm(n=100, 
            meanlog=2, 
            sdlog=0.1)

qqPlot(X, 
       distribution="lnorm", 
       meanlog=2, 
       sdlog=0.1)


###################################################
### PlotSimulated50qqPlotCar
###################################################
set.seed(123)

X <- rnorm(50, 10, 1)

library(car)

qqPlot(X, col.lines=1)

rm(X)


###################################################
### PlotSimulated100qqPlotCarLogNormal
###################################################
set.seed(123)

X <- rlnorm(100)

qqPlot(X, col.lines=1, distribution="lnorm")

rm(X)


###################################################
### PlotSimulated100hist2Normal
###################################################
set.seed(123)

X <- rnorm(50, 10)

X <- c(X, rnorm(50, 15))

hist(X, xlab=expression(x), main="")

rm(X)


###################################################
### PlotSimulated100qqplot2Normal
###################################################
set.seed(123)

X <- rnorm(50, 10)

X <- c(X, rnorm(50, 15))

qqPlot(X, col.lines=1)

rm(X)


###################################################
### Chap04.Rnw:1733-1734
###################################################
options(warn=-1)


###################################################
### KsTest
###################################################
data(OTURB)                        #
                                   #
ks.test(x=OTURB,                   # 
        y="pnorm", 
        mean = mean(OTURB), 
        sd = sd(OTURB), 
        alternative="two.sided")


###################################################
### Chap04.Rnw:1747-1748
###################################################
options(warn=0)


###################################################
### PlotBeta8020Distribution
###################################################
X <- seq(0, 1, length.out=200)

plot(X, dbeta(X, shape1=80, shape2=20), 
     type="l", 
     xlab=expression(x),
     ylab=expression("p.d.f."))

rm(X)


###################################################
### PlotBeta82Distribution
###################################################
X <- seq(0, 1, length.out=200)

plot(X, dbeta(X, shape1=8, shape2=2), 
     type="l", 
     xlab=expression(x),
     ylab=expression("p.d.f."))

rm(X)


###################################################
### PlotBeta82PosteriorDistribution
###################################################
X <- seq(0, 1, length.out=200)

plot(X, dbeta(X, shape1=8+8, shape2=2+2), 
     type="l", 
     xlab=expression(x),
     ylab=expression("p.d.f."))

lines(X, dbeta(X, shape1=8+7, shape2=2+3), 
      type="l", lty="dashed")

lines(X, dbeta(X, shape1=8+6, shape2=2+4), 
      type="l", lty="dotdash")

arrows(x1=c(0.8, 0.7, 0.6), 
       y1=c(
         dbeta(0.8, shape1=8+8, shape2=2+2),
         dbeta(0.7, shape1=8+7, shape2=2+3),
         dbeta(0.6, shape1=8+6, shape2=2+4)), 
       x0=c(0.4, 0.4, 0.4),
       y0=c(
         dbeta(0.8, shape1=8+8, shape2=2+2),
         dbeta(0.7, shape1=8+7, shape2=2+3),
         dbeta(0.6, shape1=8+6, shape2=2+4)))

text(c(0.3, 0.3, 0.3), 
     c(
       dbeta(0.8, shape1=8+8, shape2=2+2),
       dbeta(0.7, shape1=8+7, shape2=2+3),
       dbeta(0.6, shape1=8+6, shape2=2+4)), 
     c(expression(X==8), 
       expression(X==7), 
       expression(X==6)))

rm(X)


###################################################
### PlotBayesRiskFunction
###################################################
plot(c(0, 1), c(0, 5), type="n", xlab=expression(pi), ylab="Risk", ylim=c(-0.3, 5.3))

lines(c(0, 1), c(0, 5), lty="longdash")

lines(c(0, 1), c(1, 0), lty="dashed")

lines(rep(1/6, 2), c(0, 5*1/6), lty="dashed")

lines(c(0, 1/6), c(0, 5*1/6), lwd=2)

lines(c(1/6, 1), c(5*1/6, 0), lwd=2)

text(-0.01, 1.2, labels=expression(r[0]))

text(1.01, 4.8, labels=expression(r[1]))

text(1/6, -0.18, labels=expression(pi^"*"))


###################################################
### PlotPosteriorPdfAndCredibilityIntervals
###################################################
set.seed(123)
Pois <- rpois(3, lambda=2)
X <- seq(0, 6, length.out=600)

V0 <- 3
B0 <- 1
#plot(X, dgamma(X, shape=V0, scale=B0), type="l")


V1 <- V0 + sum(Pois)
B1 <- B0/(1+length(Pois)*B0)


plot(X, dgamma(X, shape=V1, scale=B1), type="l", 
     xlab=expression(x), ylab=expression(h))

ConfX <- B1/2*qchisq(p=c(0.1, 0.9), df=2*V1)
ConfY <- dgamma(ConfX, shape=V1, scale=B1)

polygon(c(X[X <= ConfX[1]], ConfX[1]), 
        c(dgamma(X[X <= ConfX[1]], shape=V1, scale=B1), 0), 
        col="grey80")
polygon(c(X[X >= ConfX[2]], ConfX[2]), 
        c(dgamma(X[X >= ConfX[2]], shape=V1, scale=B1), 0), 
        col="grey80")

rm(X, Pois, B0, B1, V0, V1, ConfX, ConfY)


###################################################
### BootstrapStudentizedTest
###################################################
library(boot)                          

data(HYBRID1)                          

set.seed(123)                          

boot.ci(boot(data=HYBRID1,             
             statistic=function(x, i) 
               mean(x[i]), 
             R=500), 
        type = "perc")                 

t.test(HYBRID1, mu=2150)               

set.seed(123)                          

B <- boot(data=HYBRID1,                
          statistic=function(x, i, mu) 
            t.test(x[i], 
                   mu=mu)$p.value,
          R=500,
          mu=2150)                     

sum(B$t <                              
      t.test(HYBRID1, 
             mu=2150)$p.value) /
  nrow(B$t)


###################################################
### BootstrapStudentizedTest2Sample
###################################################
data(HYBRID2)                          

t.test(HYBRID2$hyb1, HYBRID2$hyb2)     

set.seed(123)                          

boot(data=HYBRID2,                     
     statistic=function(x, i)               
       t.test(x=x[i,1],                
              y=x[i,2])$p.value,       
     R=500)                            


###################################################
### PlotBootstrapStudentizedTest2Sample
###################################################
data(HYBRID2)                          

set.seed(123)                          

B <- boot(data=HYBRID2,                
          statistic=function(x, i) 
            t.test(x=x[i,1], 
                   y=x[i,2])$statistic,
          R=500)                       

hist(B$t,                              
     xlab="Studentized differences", 
     main="")


###################################################
### Vartest
###################################################
data(HYBRID)

set.seed(123)

B <- apply(HYBRID, MARGIN=2, 
           FUN=boot, 
           statistic=function(x, i){
             var(x[i])
           }, 
           R = 500)

Bt0 <- sapply(B, 
              FUN=function(x) x$t0)

Bt <-  sapply(B, 
              FUN=function(x) x$t)

Bf <- max(Bt0)/min(Bt0)

FBoot <- apply(Bt, MARGIN=1, 
               FUN=function(x){
                 max(x)/min(x)
               })

Bf

quantile(FBoot, 0.95)

sum(FBoot >= Bf)/length(FBoot)

rm(Bt0, Bt, Bf, FBoot)


###################################################
### Anovatest
###################################################
onewayTestBoot <- function(x, i){
  x <- x[i,]
  y <- stack(x)
  names(y) <- c("v", "g")
  oneway.test(v ~ g, 
              data=y, 
              var.equal=TRUE)$statistic
}

set.seed(123)                          

B <- boot(data=HYBRID, 
          statistic=onewayTestBoot, 
          R=500)

B$t0

sum(B$t > B$t0)/nrow(B$t)


###################################################
### PlotAnovatest
###################################################
set.seed(123)                          

B <- boot(data=HYBRID, 
          statistic=onewayTestBoot, 
          R=500)

hist(B$t, xlab="F* values", main="")


###################################################
### BernulliSample
###################################################
qbinom(p=c(0.025, 0.975), size=50, prob=0.1)


###################################################
### BootBernulliToleranceInterval
###################################################
data(OELECT)                        

ELECINDX <- ifelse(                 
  test=OELECT >= 216 & 
    OELECT <= 224, 
  yes=1, no=0)                      

qbinomBoot <- function(x, i,        
                       size,        
                       probs=c(0.025, 
                               0.975)){
  qbinom(p=probs, 
         size=size, 
         prob=mean(x[i]))
}                                   

set.seed(123)                       

B <- boot(data=ELECINDX,            
          statistic=qbinomBoot,  
          R = 500, size = 100)   

quantile(x=B$t[,1],                  
         probs=c(0.025, 0.975))     


###################################################
### BootContinuousToleranceInterval
###################################################
data(CYCLT)

set.seed(123)                                 

B <- boot(CYCLT, 
          statistic=function(x, i){
            quantile(x[i], 
                     probs=c(0.025, 0.975))}, 
          R=500)

quantile(x=B$t[,1], probs=0.025)

quantile(x=B$t[,2], probs=0.975)


###################################################
### RandTestTwoMean
###################################################
data(OELECT1)

randomizationTest(list(a=OELECT, b=OELECT1), 
                  R=500, calc=mean, 
                  fun=function(x) x[1]-x[2],
                  seed=123)


###################################################
### Wilcox
###################################################
X <- c(0.188, 0.353, -0.257, 0.220, 0.168)

Y <- c(1.240, 1.821, 2.500, 2.319, 2.190)

wilcox.test(x=X, y=Y,
            conf.int = TRUE)

rm(X, Y)


###################################################
### Chap04End
###################################################
rm(HYBRID, HYBRID1, HYBRID2, B, CYCLT, ELECINDX, OELECT, 
   OELECT1, OTURB, onewayTestBoot, qbinomBoot)
detach(package:boot)
detach(package:AcceptanceSampling)
detach(package:car)
detach(package:mistat)
