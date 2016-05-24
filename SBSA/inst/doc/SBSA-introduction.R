### R code from vignette source 'SBSA-introduction.Rnw'

###################################################
### code chunk number 1: SBSA-introduction.Rnw:31-37
###################################################
options(width=60)
prettyVersion <- packageDescription("SBSA")$Version
prettyDate <- format(Sys.Date(), "%B %e, %Y")
library(xtable)
library(MASS)
library(SBSA)


###################################################
### code chunk number 2: SBSA-introduction.Rnw:150-156
###################################################
  set.seed(42)
  n <- 100
  tmp <- sqrt(0.6) * matrix(rnorm(n), n, 4) +
    sqrt(1 - 0.6) * matrix(rnorm(n * 4), n, 4)
  x <- tmp[, 1]
  z <- tmp[, 2:4]


###################################################
### code chunk number 3: SBSA-introduction.Rnw:160-161
###################################################
  y <- rnorm(n, x + z%*%rep(.5,3), .5)


###################################################
### code chunk number 4: SBSA-introduction.Rnw:166-169
###################################################
  w <- z[, 1:2]
  w[, 1] <- w[, 1] + rnorm(n, sd = sqrt(1/0.7 - 1))
  w[, 2] <- w[, 2] + rnorm(n, sd = sqrt(1/0.7 - 1))


###################################################
### code chunk number 5: SBSA-introduction.Rnw:173-176
###################################################
  standardize <- function(x) (x-mean(x))/sqrt(var(x))
  x.sdz <- standardize(x)
  w.sdz <- apply(w, 2, standardize)


###################################################
### code chunk number 6: SBSA-introduction.Rnw:203-205
###################################################
  a <- 6
  b <- 21


###################################################
### code chunk number 7: SBSA-introduction.Rnw:212-213
###################################################
1 - (a-1)/(a+b-2)


###################################################
### code chunk number 8: SBSA-introduction.Rnw:219-220
###################################################
pbeta(0.6, a, b)


###################################################
### code chunk number 9: SBSA-introduction.Rnw:240-243
###################################################
 sampler.jump <- c(alpha=.1, beta.z=.1,
                   sigma.sq=.1, tau.sq=.1,
                   beta.u.gamma.x=.1, gamma.z=.1)


###################################################
### code chunk number 10: SBSA-introduction.Rnw:252-254
###################################################
sbsa.fit <- fitSBSA(y, x.sdz, w.sdz, a, b, nrep = 20000,
                    sampler.jump = sampler.jump)


###################################################
### code chunk number 11: SBSA-introduction.Rnw:269-270
###################################################
names(sbsa.fit)


###################################################
### code chunk number 12: SBSA-introduction.Rnw:282-283
###################################################
sbsa.fit$acc


###################################################
### code chunk number 13: SBSA-introduction.Rnw:312-329
###################################################
sbsa.fit <- fitSBSA(y, x.sdz, w.sdz, a, b, nrep=20000, 
                    sampler.jump=c(alpha=.2, beta.z=.1,
                      sigma.sq=.1, tau.sq=.1,
                      beta.u.gamma.x=.1, gamma.z=.1))
sbsa.fit$acc

sbsa.fit <- fitSBSA(y, x.sdz, w.sdz, a, b, nrep=20000, 
                    sampler.jump=c(alpha=.15, beta.z=.1,
                      sigma.sq=.1, tau.sq=.1,
                      beta.u.gamma.x=.1, gamma.z=.1))
sbsa.fit$acc

sbsa.fit <- fitSBSA(y, x.sdz, w.sdz, a, b, nrep=20000, 
                    sampler.jump=c(alpha=.15, beta.z=.2,
                      sigma.sq=.1, tau.sq=.1,
                      beta.u.gamma.x=.1, gamma.z=.1))
sbsa.fit$acc


###################################################
### code chunk number 14: SBSA-introduction.Rnw:333-338
###################################################
sbsa.fit <- fitSBSA(y, x.sdz, w.sdz, a, b, nrep=20000, 
                    sampler.jump=c(alpha=.15, beta.z=.2,
                      sigma.sq=.35, tau.sq=.1,
                      beta.u.gamma.x=.7, gamma.z=1.1))
sbsa.fit$acc


###################################################
### code chunk number 15: traceplotCode
###################################################
mfrow <- par(mfrow=c(2,2))
plot(window(ts(sbsa.fit$alpha[,1]), deltat=30), ylab=expression(alpha[0]))
plot(window(ts(sbsa.fit$alpha[,2]), deltat=30), ylab=expression(alpha[x]))
plot(window(ts(sbsa.fit$beta.u), deltat=30), ylab=expression(beta[u]))
plot(window(ts(sbsa.fit$gamma.x), deltat=30), ylab=expression(gamma[x]))
par(mfrow=mfrow)


###################################################
### code chunk number 16: traceplotFig
###################################################
mfrow <- par(mfrow=c(2,2))
plot(window(ts(sbsa.fit$alpha[,1]), deltat=30), ylab=expression(alpha[0]))
plot(window(ts(sbsa.fit$alpha[,2]), deltat=30), ylab=expression(alpha[x]))
plot(window(ts(sbsa.fit$beta.u), deltat=30), ylab=expression(beta[u]))
plot(window(ts(sbsa.fit$gamma.x), deltat=30), ylab=expression(gamma[x]))
par(mfrow=mfrow)


###################################################
### code chunk number 17: SBSA-introduction.Rnw:375-377
###################################################
  mean(sbsa.fit$alpha[10001:20000, 2])
  sqrt(var(sbsa.fit$alpha[10001:20000, 2]))


###################################################
### code chunk number 18: SBSA-introduction.Rnw:384-386
###################################################
trgt <- sbsa.fit$alpha[10001:20000,2]/sqrt(var(x))
c(mean(trgt),sqrt(var(trgt)))


###################################################
### code chunk number 19: SBSA-introduction.Rnw:401-406
###################################################
  w <- z[, 1:2]
  w[, 1] <- w[, 1] + rnorm(n, sd = sqrt(1/0.7 - 1))
  w[, 2] <- w[, 2] + rnorm(n, sd = sqrt(1/0.95 - 1))

  w.sdz <- apply(w, 2, standardize)


###################################################
### code chunk number 20: SBSA-introduction.Rnw:416-418
###################################################
  a <- c(6, 3)
  b <- c(21, 39)


###################################################
### code chunk number 21: SBSA-introduction.Rnw:423-424
###################################################
1 - (a-1)/(a+b-2)


###################################################
### code chunk number 22: SBSA-introduction.Rnw:428-429
###################################################
pbeta(c(0.6, 0.8), a, b)


###################################################
### code chunk number 23: SBSA-introduction.Rnw:443-448
###################################################
sbsa.fit <- fitSBSA(y, x.sdz, w.sdz, a, b, nrep=20000, 
                    sampler.jump=list(alpha=.15, beta.z=.2,
                      sigma.sq=.35, tau.sq=c(.1, .05),
                      beta.u.gamma.x=.7, gamma.z=1.1))
sbsa.fit$acc


###################################################
### code chunk number 24: SBSA-introduction.Rnw:453-458
###################################################
sbsa.fit <- fitSBSA(y, x.sdz, w.sdz, a, b, nrep=20000, 
                    sampler.jump=list(alpha=.14, beta.z=.15,
                      sigma.sq=.25, tau.sq=c(.1, .05),
                      beta.u.gamma.x=.6, gamma.z=1.1))
sbsa.fit$acc


###################################################
### code chunk number 25: traceplotCode
###################################################
mfrow <- par(mfrow=c(2,2))
plot(window(ts(sbsa.fit$alpha[,1]), deltat=30), ylab=expression(alpha[0]))
plot(window(ts(sbsa.fit$alpha[,2]), deltat=30), ylab=expression(alpha[x]))
plot(window(ts(sbsa.fit$beta.u), deltat=30), ylab=expression(beta[u]))
plot(window(ts(sbsa.fit$gamma.x), deltat=30), ylab=expression(gamma[x]))
par(mfrow=mfrow)


###################################################
### code chunk number 26: traceplotFig
###################################################
mfrow <- par(mfrow=c(2,2))
plot(window(ts(sbsa.fit$alpha[,1]), deltat=30), ylab=expression(alpha[0]))
plot(window(ts(sbsa.fit$alpha[,2]), deltat=30), ylab=expression(alpha[x]))
plot(window(ts(sbsa.fit$beta.u), deltat=30), ylab=expression(beta[u]))
plot(window(ts(sbsa.fit$gamma.x), deltat=30), ylab=expression(gamma[x]))
par(mfrow=mfrow)


###################################################
### code chunk number 27: tauCode
###################################################
tau.density <- kde2d(sbsa.fit$tau.sq[, 1], 
                     sbsa.fit$tau.sq[, 2], 
                     lims = c(0, max(sbsa.fit$tau), 
                              0, max(sbsa.fit$tau)))
filled.contour(tau.density,
               color.palette = function(n) grey(n:0 / n),
               xlab = expression({tau[1]}^2), 
               ylab = expression({tau[2]}^2))


###################################################
### code chunk number 28: tauPlot
###################################################
tau.density <- kde2d(sbsa.fit$tau.sq[, 1], 
                     sbsa.fit$tau.sq[, 2], 
                     lims = c(0, max(sbsa.fit$tau), 
                              0, max(sbsa.fit$tau)))
filled.contour(tau.density,
               color.palette = function(n) grey(n:0 / n),
               xlab = expression({tau[1]}^2), 
               ylab = expression({tau[2]}^2))


###################################################
### code chunk number 29: SBSA-introduction.Rnw:509-511
###################################################
  mean(sbsa.fit$alpha[10001:20000, 2])
  sqrt(var(sbsa.fit$alpha[10001:20000, 2]))


###################################################
### code chunk number 30: SBSA-introduction.Rnw:516-518
###################################################
trgt <- sbsa.fit$alpha[10001:20000,2]/sqrt(var(x))
c(mean(trgt),sqrt(var(trgt)))


