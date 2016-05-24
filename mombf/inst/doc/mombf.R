### R code from vignette source 'mombf.Rnw'

###################################################
### code chunk number 1: fig1plot
###################################################
 library(mombf)
 tau <- 1
 thseq <- seq(-3,3,length=1000)
 plot(thseq,dmom(thseq,tau=tau),type='l',ylab='Prior density')
 lines(thseq,dmom(thseq,tau=tau,baseDensity='t',penalty='quadratic',nu=3),lty=2,col=2)
 lines(thseq,dimom(thseq,tau=tau),lty=3,col=3)


###################################################
### code chunk number 2: fig1
###################################################
 library(mombf)
 tau <- 1
 thseq <- seq(-3,3,length=1000)
 plot(thseq,dmom(thseq,tau=tau),type='l',ylab='Prior density')
 lines(thseq,dmom(thseq,tau=tau,baseDensity='t',penalty='quadratic',nu=3),lty=2,col=2)
 lines(thseq,dimom(thseq,tau=tau),lty=3,col=3)


###################################################
### code chunk number 3: fig1bplot
###################################################
 library(mombf)
 plot(thseq,pmom(thseq,tau=tau),type='l',ylab='Prior cdf')
 lines(thseq,pimom(thseq,tau=tau),lty=3,col=3)


###################################################
### code chunk number 4: fig1b
###################################################
 library(mombf)
 plot(thseq,pmom(thseq,tau=tau),type='l',ylab='Prior cdf')
 lines(thseq,pimom(thseq,tau=tau),lty=3,col=3)


###################################################
### code chunk number 5: one
###################################################
data(hald)
dim(hald)
lm1 <- lm(hald[,1] ~ hald[,2] + hald[,3] + hald[,4] + hald[,5])
summary(lm1)


###################################################
### code chunk number 6: two
###################################################
prior.mode <- .2^2
V <- summary(lm1)$cov.unscaled
diag(V)
taumom <- mode2g(prior.mode,prior='normalMom')
tautmom <- mode2g(prior.mode,prior='tMom',nu=3)
tauimom <- mode2g(prior.mode,prior='iMom')
taumom
tautmom
tauimom


###################################################
### code chunk number 7: fig2plot
###################################################
thseq <- seq(-1,1,length=1000)
plot(thseq,dmom(thseq,V1=nrow(hald)*V[2,2],tau=taumom),type='l',xlab='theta/sigma',ylab='Prior density')
lines(thseq,dmom(thseq,V1=nrow(hald)*V[2,2],tau=tautmom,baseDensity='t',nu=3,penalty='quadratic'),lty=2,col=2)
lines(thseq,dimom(thseq,V1=nrow(hald)*V[2,2],tau=tauimom),lty=3,col=3)
abline(v=.2,lty=2,col='gray')


###################################################
### code chunk number 8: fig2
###################################################
thseq <- seq(-1,1,length=1000)
plot(thseq,dmom(thseq,V1=nrow(hald)*V[2,2],tau=taumom),type='l',xlab='theta/sigma',ylab='Prior density')
lines(thseq,dmom(thseq,V1=nrow(hald)*V[2,2],tau=tautmom,baseDensity='t',nu=3,penalty='quadratic'),lty=2,col=2)
lines(thseq,dimom(thseq,V1=nrow(hald)*V[2,2],tau=tauimom),lty=3,col=3)
abline(v=.2,lty=2,col='gray')


###################################################
### code chunk number 9: twobis
###################################################
a <- .2; priorp <- .05
taumom2 <- priorp2g(priorp=priorp,q=a,prior='normalMom')
tauimom2 <- priorp2g(priorp=priorp,q=-a,prior='iMom')
taumom2
tauimom2


###################################################
### code chunk number 10: three
###################################################
set.seed(4*2*2008)
mombf(lm1,coef=2,g=taumom)
mombf(lm1,coef=2,g=tautmom,baseDensity='t')
imombf(lm1,coef=2,g=tauimom,method='adapt')
imombf(lm1,coef=2,g=tauimom,method='MC',B=10^5)
zellnerbf(lm1,coef=2,g=1)


###################################################
### code chunk number 11: four
###################################################
imombf(lm1,coef=2,g=tauimom,method='MC',B=10^5)


###################################################
### code chunk number 12: five
###################################################
sr <- sqrt(sum(lm1$residuals^2)/(nrow(hald)-5))
thest <- coef(lm1)[2]/sr
thest


###################################################
### code chunk number 13: fig3plot
###################################################
prior.mode <- seq(.01,1,length=100)^2
taumom <- mode2g(prior.mode,prior='normalMom')
bf1 <- mombf(lm1,coef=2,g=taumom)
bf2 <- zellnerbf(lm1,coef=2,g=taumom)
plot(prior.mode,bf1,type='l',ylab='BF',ylim=range(c(bf1,bf2)))
lines(prior.mode,bf2,lty=2,col=2)
abline(v=thest,lty=2)


###################################################
### code chunk number 14: fig3
###################################################
prior.mode <- seq(.01,1,length=100)^2
taumom <- mode2g(prior.mode,prior='normalMom')
bf1 <- mombf(lm1,coef=2,g=taumom)
bf2 <- zellnerbf(lm1,coef=2,g=taumom)
plot(prior.mode,bf1,type='l',ylab='BF',ylim=range(c(bf1,bf2)))
lines(prior.mode,bf2,lty=2,col=2)
abline(v=thest,lty=2)


###################################################
### code chunk number 15: six
###################################################
set.seed(4*2*2008)
n <- 50; theta <- c(log(2),0)
x <- matrix(NA,nrow=n,ncol=2)
x[,1] <- rnorm(n,0,1); x[,2] <- rnorm(n,.5*x[,1],1)
p <- pnorm(x %*% matrix(theta,ncol=1))
y <- rbinom(n,1,p)


###################################################
### code chunk number 16: seven
###################################################
glm1 <- glm(y~x[,1]+x[,2],family=binomial(link = "probit"))
thetahat <- coef(glm1)
V <- summary(glm1)$cov.scaled


###################################################
### code chunk number 17: eight
###################################################
g <- .5
bfmom.1 <- momknown(thetahat[2],V[2,2],n=n,g=g,sigma=1)
bfimom.1 <- imomknown(thetahat[2],V[2,2],n=n,nuisance.theta=2,g=g,sigma=1)
bfmom.1
bfimom.1


###################################################
### code chunk number 18: nine
###################################################
bfmom.2 <- momknown(thetahat[3],V[3,3],n=n,g=g,sigma=1)
bfimom.2 <- imomknown(thetahat[3],V[3,3],n=n,nuisance.theta=2,g=g,sigma=1)
bfmom.2
bfimom.2


###################################################
### code chunk number 19: ten
###################################################
bfmom.0 <- momknown(thetahat[2:3],V[2:3,2:3],n=n,g=g,sigma=1)
bfimom.0 <- imomknown(thetahat[2:3],V[2:3,2:3],n=n,nuisance.theta=2,g=g,sigma=1)
bfmom.0
bfimom.0


###################################################
### code chunk number 20: eleven
###################################################
prior.prob <- rep(1/4,4)
bf <- c(bfmom.0,bfmom.1,bfmom.2,1)
pos.prob <- prior.prob*bf/sum(prior.prob*bf)
pos.prob


###################################################
### code chunk number 21: varsel1
###################################################
set.seed(2011*01*18)
x <- matrix(rnorm(100*3),nrow=100,ncol=3)
theta <- matrix(c(1,1,0),ncol=1)
y <- x %*% theta + rnorm(100)


###################################################
### code chunk number 22: varsel2
###################################################
priorCoef <- imomprior(tau=.131)
priorDelta <- modelbbprior(alpha.p=1,beta.p=1)
priorVar <- igprior(alpha=.01,lambda=.01)


###################################################
### code chunk number 23: varsel3
###################################################
fit1 <- modelSelection(y=y, x=x, center=FALSE, scale=FALSE, niter=10^2,
priorCoef=priorCoef, priorDelta=priorDelta, priorVar=priorVar, 
method='Laplace')
fit1$postMode
fit1$margpp


###################################################
### code chunk number 24: varsel4
###################################################
correct <- t(fit1$postSample)==c(TRUE,TRUE,FALSE)
table(colSums(correct)==3)


