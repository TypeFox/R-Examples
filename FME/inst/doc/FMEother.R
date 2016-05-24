### R code from vignette source 'FMEother.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library("FME")
options(prompt = "> ")
options(width=70)


###################################################
### code chunk number 2: FMEother.Rnw:110-111
###################################################
require(FME)


###################################################
### code chunk number 3: FMEother.Rnw:114-116
###################################################
Obs <- data.frame(x=c(   28,  55,   83,  110,  138,  225,  375),   # mg COD/l
                  y=c(0.053,0.06,0.112,0.105,0.099,0.122,0.125))   # 1/hour


###################################################
### code chunk number 4: FMEother.Rnw:119-120
###################################################
Model <- function(p, x) return(data.frame(x = x, y = p[1]*x/(x+p[2])))


###################################################
### code chunk number 5: FMEother.Rnw:127-128
###################################################
Residuals  <- function(p) (Obs$y - Model(p, Obs$x)$y)


###################################################
### code chunk number 6: FMEother.Rnw:132-135
###################################################
print(system.time(
P      <- modFit(f = Residuals, p = c(0.1, 1))
))


###################################################
### code chunk number 7: FMEother.Rnw:139-141
###################################################
sP    <- summary(P)
sP


###################################################
### code chunk number 8: FMEother.Rnw:145-146
###################################################
x      <-0:375


###################################################
### code chunk number 9: Monplot
###################################################
par(mfrow = c(2, 2))
plot(P, mfrow = NULL)
plot(Obs, pch = 16, cex = 2, xlim = c(0, 400), ylim = c(0, 0.15),
     xlab = "mg COD/l", ylab = "1/hr", main = "best-fit")
lines(Model(P$par, x))
par(mfrow = c(1, 1))


###################################################
### code chunk number 10: Monplot
###################################################
par(mfrow = c(2, 2))
plot(P, mfrow = NULL)
plot(Obs, pch = 16, cex = 2, xlim = c(0, 400), ylim = c(0, 0.15),
     xlab = "mg COD/l", ylab = "1/hr", main = "best-fit")
lines(Model(P$par, x))
par(mfrow = c(1, 1))


###################################################
### code chunk number 11: FMEother.Rnw:180-187
###################################################
Covar   <- sP$cov.scaled * 2.4^2/2
s2prior <- sP$modVariance

print(system.time(
MCMC <- modMCMC(f = Residuals, p = P$par, jump = Covar, niter = 3000,
                var0 = s2prior, wvar0 = 1, lower = c(0, 0))
))


###################################################
### code chunk number 12: FMEother.Rnw:191-196
###################################################
print(system.time(
MCMC <- modMCMC(f = Residuals, p = P$par, jump = Covar, niter = 3000, 
      ntrydr = 3, var0 = s2prior, wvar0 = 1, updatecov = 100, lower = c(0, 0))
))
MCMC$count


###################################################
### code chunk number 13: Monmcmc
###################################################
plot(MCMC, Full = TRUE)


###################################################
### code chunk number 14: Monmcmc
###################################################
plot(MCMC, Full = TRUE)


###################################################
### code chunk number 15: Monhist
###################################################
hist(MCMC, Full = TRUE, col = "darkblue")


###################################################
### code chunk number 16: Monhist
###################################################
hist(MCMC, Full = TRUE, col = "darkblue")


###################################################
### code chunk number 17: Monpairs
###################################################
pairs(MCMC)


###################################################
### code chunk number 18: Monpairs
###################################################
pairs(MCMC)


###################################################
### code chunk number 19: FMEother.Rnw:244-247
###################################################
cor(MCMC$pars)
cov(MCMC$pars)
sP$cov.scaled


###################################################
### code chunk number 20: FMEother.Rnw:253-255
###################################################
MC <- as.mcmc(MCMC$pars)
raftery.diag(MC)


###################################################
### code chunk number 21: cumu
###################################################
cumuplot(MC)


###################################################
### code chunk number 22: cumu
###################################################
cumuplot(MC)


###################################################
### code chunk number 23: FMEother.Rnw:278-279
###################################################
sR<-sensRange(parInput=MCMC$pars,func=Model,x=1:375)


###################################################
### code chunk number 24: Monsum
###################################################
plot(summary(sR), quant = TRUE)
points(Obs)


###################################################
### code chunk number 25: Monsum
###################################################
plot(summary(sR), quant = TRUE)
points(Obs)


###################################################
### code chunk number 26: FMEother.Rnw:312-313
###################################################
pset <- attributes(sR)$pset


###################################################
### code chunk number 27: FMEother.Rnw:318-323
###################################################
nout  <- nrow(sR)
sR2   <- sR
ivar  <- 3:ncol(sR)
error <- rnorm(nout, mean = 0, sd = sqrt(MCMC$sig[pset]))
sR2[,ivar] <- sR2[ ,ivar] + error


###################################################
### code chunk number 28: MonsumM
###################################################
plot(summary(sR2),quant=TRUE)
points(Obs)


###################################################
### code chunk number 29: MonsumM
###################################################
plot(summary(sR2),quant=TRUE)
points(Obs)


