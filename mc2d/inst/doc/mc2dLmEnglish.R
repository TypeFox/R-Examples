### R code from vignette source 'mc2dLmEnglish.rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: sh0
###################################################
options("width"=100,"digits"=3)
set.seed(666)


###################################################
### code chunk number 2: sh1
###################################################
Nmax <- 7.3
murefLm <- 6.2; TminLm <- -2.9

murefFF <- 4.1; TminFF <- -4.5

Lm0 <- log10(1); FF0 <- 2.78

d1 <- 1.1; mT1 <- 3.2; sdT1 <- 2.1

d2 <- 4.7; mT2 <- 5.5; sdT2 <- 1.0

d3 <- 4.3; mT3 <- 8.2; sdT3 <- 2.0

conso <- 35
r <- 4.7e-14


modGrowth <- function(duration, mTemp, sdTemp,
                      N0Lm, murefLm, TminLm,
                      N0FF, murefFF, TminFF,
                      Nmax, Tref=25) {
  N0Lm <- pmin(N0Lm, Nmax)
  N0FF <- pmin(N0FF, Nmax)
  dLm <- (Nmax-N0Lm) * log(10)/murefLm * (Tref-TminLm)^2 / (sdTemp^2 + (mTemp-TminLm)^2)
  dLm <- ifelse(mTemp < TminLm & N0Lm!=Nmax, Inf, dLm)
  dFF <- (Nmax-N0FF) * log(10)/murefFF * (Tref-TminFF)^2 / (sdTemp^2 + (mTemp-TminFF)^2)
  dFF <- ifelse(mTemp < TminFF & N0FF!=Nmax, Inf, dFF)
  realDuration <- pmin(duration, dLm , dFF)
  xLm <- N0Lm + (mTemp > TminLm) * murefLm/log(10) * 
        (sdTemp^2 + (mTemp - TminLm)^2) /  ((Tref - TminLm)^2) * realDuration
  xFF <- N0FF + (mTemp > TminFF) * murefFF/log(10) * 
        (sdTemp^2 + (mTemp - TminFF)^2) /  ((Tref - TminFF)^2) * realDuration
  return(list(xLm = xLm, xFF=xFF))}

x1 <- modGrowth(d1, mT1, sdT1,
                Lm0, murefLm, TminLm,
                FF0, murefFF, TminFF,
                Nmax)
x2 <- modGrowth(d2, mT2, sdT2,
                x1$xLm, murefLm, TminLm,
                x1$xFF, murefFF, TminFF,
                Nmax)
x3 <- modGrowth(d3, mT3, sdT3,
                x2$xLm, murefLm, TminLm,
                x2$xFF, murefFF, TminFF,
                Nmax)
x3
conta <-10^x3$xLm
conta
expo <- conso * conta
expo
risk <- 1 - (1 - r)^expo
risk


###################################################
### code chunk number 3: shCallLibrary
###################################################
library(fitdistrplus)
library(mc2d)
ndvar(10001)


###################################################
### code chunk number 4: shLm0FF0V
###################################################
dataC <- data.frame(
  left =c(rep(NA,43), rep(.2,7),.3,rep(.4,4),1,1.6,.6,.6,2.4,5.4,7),
  right=c(rep(0.2,43),rep(.2,7),.3,rep(.4,4),1,1.6,.6,.6,2.4,5.4,7)
  )
fit <- fitdistcens(log10(dataC), "norm")
fit
Lm0V <- mcstoc(rnorm, mean = fit$est["mean"], sd = fit$est["sd"], rtrunc=TRUE, linf=-2)
FF0V <- mcstoc(rnorm, mean=2.78, sd=1.14)


###################################################
### code chunk number 5: shGrowtV
###################################################

NmaxV    <- mcstoc(rnorm, mean=7.27, sd = 0.86)

murefLmV <- mcstoc(rnorm, mean = 6.24, sd = 0.75, rtrunc=TRUE, linf=0)
TminLmV  <- mcstoc(rnorm, mean = -2.86, sd = 1.93)

murefFFV <- mcstoc(rnorm, mean = 4.12, sd = 1.97, rtrunc=TRUE, linf=0)
TminFFV <-  mcstoc(rnorm, mean = -4.52, sd = 7.66)


###################################################
### code chunk number 6: shTtV
###################################################
d1V <- mcstoc(rexp, rate = 1/1.1)
mT1V <- mcstoc(rnorm, mean = 3.2, sd = 2.2, rtrunc = TRUE, linf = -3, lsup = 25)
sdT1V <- sqrt(mcstoc(rgamma, shape = 1.16, scale=4.61))

d2V <- mcstoc(rexp, rate = 1/4.7, rtrunc=TRUE, lsup=28-d1V)
mT2V <- mcstoc(rnorm, mean = 5.5, sd = 2.2, rtrunc = TRUE, linf = -3, lsup = 25)
sdT2V <- sqrt(mcstoc(rgamma, shape = 0.65, scale=2.09))

d3V <- mcstoc(rexp, rate = 1/4.3, rtrunc=TRUE, lsup=28-(d1V+d2V))
mT3V <- mcstoc(rnorm, mean = 8.2, sd = 3.8, rtrunc = TRUE, linf = -3, lsup = 25)
sdT3V <- sqrt(mcstoc(rgamma, shape = 0.35, scale=19.7))


###################################################
### code chunk number 7: shConsoV
###################################################
consoV <- mcstoc(rempiricalD, 
	values = c(10, 12, 19, 20, 30, 34, 40, 50, 60, 67.5, 80, 100, 250), 
	prob = c(11, 1, 1, 29, 12, 1, 41, 4, 4, 1, 4, 1, 1))


###################################################
### code chunk number 8: shModelV
###################################################
r <- mcdata(4.7e-14, type = "0")
x1V <- modGrowth(d1V, mT1V, sdT1V,
                Lm0V, murefLmV, TminLmV,
                FF0V, murefFFV, TminFFV,
                NmaxV)
x2V <- modGrowth(d2V, mT2V, sdT2V,
                x1V$xLm, murefLmV, TminLmV,
                x1V$xFF, murefFFV, TminFFV,
                NmaxV)
x3V <- modGrowth(d3V, mT3V, sdT3V,
                x2V$xLm, murefLmV, TminLmV,
                x2V$xFF, murefFFV, TminFFV,
                NmaxV)

contaV <-10^x3V$xLm
expoV <- consoV * contaV
riskV <- 1 - exp(-r * expoV )
Lm1 <- mc(Lm0V, FF0V, NmaxV, murefLmV, TminLmV, murefFFV, TminFFV,
          d1V, mT1V, sdT1V, d2V, mT2V, sdT2V, d3V, mT3V, sdT3V,
          consoV, r, contaV, expoV, riskV)
Lm1
sLm1 <- mc(contaV=Lm1$contaV, expoV=Lm1$expoV, riskV=Lm1$riskV)
summary(sLm1, probs = c(0, 0.5, 0.75, 0.95, 1))


###################################################
### code chunk number 9: sh4
###################################################
meanRisk <-  mcapply(riskV,"var",mean)
expectedN <- round(0.065 * unmc(meanRisk) * 6.4 * 49090000)
expectedN


###################################################
### code chunk number 10: sh8
###################################################
ndunc(101)
bootLm0 <- bootdistcens(fit, niter=ndunc())
MLm0 <- mcdata(bootLm0$est$mean,type="U")
SLm0 <- mcdata(bootLm0$est$sd,type="U")
Lm0VU <- mcstoc(rnorm, type="VU", mean=MLm0, sd=SLm0, rtrunc=TRUE, linf=-2)


###################################################
### code chunk number 11: shFF0VU
###################################################
MLm0FF <- mcstoc(rnorm, type="U", mean=2.78, sd=0.265)
SLm0FF <- mcstoc(rlnorm, type="U", meanlog=0.114, sdlog=0.172)
FF0VU <- mcstoc(rnorm, type="VU", mean=MLm0FF, sd=SLm0FF)


###################################################
### code chunk number 12: sh5
###################################################
MmurefLm <- mcstoc(rgamma, type="U", shape=69.7, scale=0.0896)
SmurefLm <- mcstoc(rlnorm, type="U", meanlog = 1.03, sdlog = 0.191)
murefLmVU <- mcstoc(rnorm, type="VU", mean=MmurefLm, sd=SmurefLm, rtrunc=TRUE, linf=0)

MTminLm <- mcstoc(rnorm, type="U", mean=-2.86, sd=0.459)
STminLm <- mcstoc(rlnorm, type="U", meanlog = 0.638, sdlog = 0.208)
TminLmVU <- mcstoc(rnorm, type="VU", mean = MTminLm, sd = STminLm, rtrunc=TRUE, lsup=25)

MmurefFF <- mcstoc(rgamma, type="U", shape=32.5, scale=.127)
SmurefFF <- mcstoc(rlnorm, type="U", meanlog = -.656, sdlog = 0.221)
murefFFVU <- mcstoc(rnorm, type="VU", mean=MmurefFF, sd=SmurefFF, rtrunc=TRUE, linf=0)

MTminFF <- mcstoc(rnorm, type="U", mean=-4.52, sd=1.23)
STminFF <- mcstoc(rlnorm, type="U", meanlog = 2.00, sdlog = 0.257)
TminFFVU <- mcstoc(rnorm, type="VU", mean = MTminFF, sd = STminFF, rtrunc=TRUE, lsup=25)

MNmax <- mcstoc(rnorm, type="U", mean=7.27, sd=0.276)
SNmax <- mcstoc(rlnorm, type="U", meanlog = -0.172, sdlog = 0.218)
NmaxVU <- mcstoc(rnorm, type="VU", mean = MNmax, sd = SNmax)


###################################################
### code chunk number 13: shPrevU
###################################################
prevU <- mcstoc(rbeta,type="U", shape1=41+1, shape2=626-41+1)


###################################################
### code chunk number 14: sh9
###################################################
x1VU <- modGrowth(d1V, mT1V, sdT1V,
                Lm0VU, murefLmVU, TminLmVU,
                FF0VU, murefFFVU, TminFFVU,
                NmaxVU)
x2VU <- modGrowth(d2V, mT2V, sdT2V,
                x1VU$xLm, murefLmVU, TminLmVU,
                x1VU$xFF, murefFFVU, TminFFVU,
                NmaxVU)
x3VU <- modGrowth(d3V, mT3V, sdT3V,
                x2VU$xLm, murefLmVU, TminLmVU,
                x2VU$xFF, murefFFVU, TminFFVU,
                NmaxVU)

contaVU <-10^x3VU$xLm
expoVU <- consoV * contaVU
riskVU <- 1 - exp(-r * expoVU)
Lm2 <- mc(Lm0VU, FF0VU, NmaxVU, murefLmVU, TminLmVU, murefFFVU, TminFFVU,
          d1V, mT1V, sdT1V, d2V, mT2V, sdT2V, d3V, mT3V, sdT3V,
          consoV, r, contaVU, expoVU, riskVU)
Lm2
sLm2 <- mc(contaVU=Lm2$contaVU, expoVU=Lm2$expoVU, riskVU=Lm2$riskVU)
summary(sLm2, probs = c(0, 0.5, 0.75, 0.95, 1))


###################################################
### code chunk number 15: sh7
###################################################
meanRiskU <-  mcapply(riskVU,"var",mean)
expectedNU <- round(prevU * meanRiskU * 6.4 * 49090000)
summary(expectedNU)


###################################################
### code chunk number 16: sh3
###################################################
torn <- tornado(Lm2)
torn
tornunc <- tornadounc(Lm2, quant=.975)
tornunc


###################################################
### code chunk number 17: sh3b (eval = FALSE)
###################################################
## plot(torn)
## plot(tornunc, stat="mean risk")


###################################################
### code chunk number 18: sh3c
###################################################
plot(torn)


###################################################
### code chunk number 19: sh3d
###################################################
plot(tornunc, stat="mean risk")


