### R code from vignette source 'TBSSurvival.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: TBSSurvival.Rnw:108-110
###################################################
options(width=60)
set.seed(1)


###################################################
### code chunk number 2: example1
###################################################
time <- pmin(rgamma(30,10,2),rep(6,30))
delta <- rep(1,30)
for (i in 1:30) {
  if (time[i] == 6) delta[i] <- 0
}
data <- cbind(time,delta)
data


###################################################
### code chunk number 3: tbsinst (eval = FALSE)
###################################################
## install.packages("TBSSurvival_VERSION.tar.gz",
##                  repos=NULL,type="source") ## from local file
## install.packages("TBSSurvival")  ## or from CRAN


###################################################
### code chunk number 4: <tbs
###################################################
library("TBSSurvival")


###################################################
### code chunk number 5: fig1a
###################################################
x <- seq(0.1,10,0.1)
plot(x, (1-ptbs(x,lambda=1,xi=sqrt(2),beta=1,
                dist=dist.error("norm"))),
     type="l",lwd=2,lty=1,col=2,ylim=c(0,1),xlab="t",
     ylab="R(t)",main="Reliability",cex.lab=1.2)


###################################################
### code chunk number 6: fig1b
###################################################
plot(x, (dtbs(x,lambda=1,xi=sqrt(2),beta=1,
              dist=dist.error("norm"))),
     type="l",lwd=2,lty=1,col=2,xlab="t",
     ylab="f(t)",main="Density",cex.lab=1.2)


###################################################
### code chunk number 7: fig1c
###################################################
plot(x, (htbs(x,lambda=1,xi=sqrt(2),beta=1,
              dist=dist.error("norm"))),
     type="l",lwd=2,lty=1,col=2,xlab="t",
     ylab="h(t)",main="Hazard",cex.lab=1.2)


###################################################
### code chunk number 8: fig2a
###################################################
x <- seq(0.1,10,0.1)
plot(x, (1-ptbs(x, lambda=1, xi=1, beta=1, dist=dist.error("doubexp"))),
     type="l", lwd=2, lty=1, col=2, ylim=c(0,1), xlab="t", 
     ylab="R(t)", main="Reliability", cex.lab=1.2)


###################################################
### code chunk number 9: fig2b
###################################################
plot(x, (dtbs(x, lambda=1, xi=1, beta=1, dist=dist.error("doubexp"))),
     type="l", lwd=2, lty=1, col=2, xlab="t", ylab="f(t)",
     main="Density", cex.lab=1.2)


###################################################
### code chunk number 10: fig2c
###################################################
plot(x, (htbs(x, lambda=1, xi=1, beta=1, 
                    dist=dist.error("doubexp"))), type="l", lwd=2, 
     lty=1, col=2, xlab="t", ylab="h(t)",
     main="Hazard", cex.lab=1.2)


###################################################
### code chunk number 11: <mynormal
###################################################
mynormal = list(
  d = function(x,xi) dnorm(x,mean=0,sd=sqrt(xi)), # density
  p = function(x,xi) pnorm(x,mean=0,sd=sqrt(xi)), # distr
  q = function(x,xi) qnorm(x,mean=0,sd=sqrt(xi)), # quantile
  r = function(x,xi) rnorm(x,mean=0,sd=sqrt(xi)), # generation
  name = "norm"
  )


###################################################
### code chunk number 12: mle
###################################################
formula <- survival::Surv(data[,1],data[,2]) ~ 1
tbs.mle <- tbs.survreg.mle(formula,dist=mynormal, nstart=3,
                           method="Nelder-Mead")
tbs.mle


###################################################
### code chunk number 13: mlekm
###################################################
# Kaplan-Meier estimation
km <- survival::survfit(formula)
plot(tbs.mle,lwd=2,col="gray20",ylab="R(t)",
     xlab="t: number of cycles (in thousands)",
     main="Reliability function (MLE)",cex.lab=1.2)
lines(km)


###################################################
### code chunk number 14: covars (eval = FALSE)
###################################################
## library(survival)
## data(colon)
## ## Running MLE on colon (from survival package) with a covariate
## colon$age60=as.numeric(colon$age>60) #threshold defined from medical papers
## s=tbs.survreg.mle(Surv(colon$time,colon$status==1) ~ colon$age60,
##       dist=mynormal,method=c("Nelder-Mead"),nstart=3,verbose=FALSE)
## summary(s)


###################################################
### code chunk number 15: be
###################################################
tbs.be <- tbs.survreg.be(formula,dist=mynormal,
                         guess.lambda=2,guess.xi=4,
                         guess.beta=1.5,burn=1000,
                         jump=10,size=1000,scale=0.06)


###################################################
### code chunk number 16: bekm
###################################################
plot(tbs.be,ylab="R(t)",
     xlab="t: number of cycles (in thousands)",
     main="Reliability function (BE)",
     lwd=2,lty=1,col=2,lwd.HPD=2,lty.HPD=2,col.HPD=2)
lines(km)


