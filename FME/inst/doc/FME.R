### R code from vignette source 'FME.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library("FME")
options("prompt" = "R> ", "continue" = "+  ")
options(width=80)
set.seed(1257)


###################################################
### code chunk number 2: FME.Rnw:349-370
###################################################
HIV_R <- function (pars, V_0 = 50000, dV_0 = -200750, T_0 = 100) {

  derivs <- function(time, y, pars) {
    with (as.list(c(pars, y)), {
      dT <- lam - rho * T - bet * T * V
      dI <- bet * T * V - delt * I
      dV <- n * delt * I - c * V - bet * T * V

      return(list(c(dT, dI, dV), logV = log(V)))
    })
  }

  # initial conditions
  I_0   <- with(as.list(pars), (dV_0 + c * V_0) / (n * delt))
  y     <- c(T = T_0, I = I_0, V = V_0)

  times <- c(seq(0, 0.8, 0.1), seq(2, 60, 2))
  out   <- ode(y = y, parms = pars, times = times, func = derivs)

  as.data.frame(out)
}


###################################################
### code chunk number 3: FME.Rnw:384-395
###################################################
HIV <- function (pars, V_0 = 50000, dV_0 = -200750, T_0 = 100) {

  I_0 <- with(as.list(pars), (dV_0 + c * V_0) / (n * delt))
  y <- c(T = T_0, I = I_0, V = V_0)

  times <- c(0, 0.1, 0.2, 0.4, 0.6, 0.8, seq(2, 60, by = 2))
  out <- ode(y = y, parms = pars, times = times, func = "derivshiv",
    initfunc = "inithiv", nout = 1, outnames = "logV", dllname = "FME")

  as.data.frame(out)
}


###################################################
### code chunk number 4: FME.Rnw:404-406
###################################################
pars <- c(bet = 0.00002, rho = 0.15, delt = 0.55, c = 5.5, lam = 80, n = 900)
out <- HIV(pars = pars)


###################################################
### code chunk number 5: ode
###################################################
par(mfrow = c(1, 2))
plot(out$time, out$logV, main = "Viral load", ylab = "log(V)",
  xlab = "time", type = "b")
plot(out$time, out$T, main = "CD4+ T", ylab = "-", xlab = "time", type = "b")
par(mfrow = c(1, 1))


###################################################
### code chunk number 6: odefig
###################################################
par(mfrow = c(1, 2))
plot(out$time, out$logV, main = "Viral load", ylab = "log(V)",
  xlab = "time", type = "b")
plot(out$time, out$T, main = "CD4+ T", ylab = "-", xlab = "time", type = "b")
par(mfrow = c(1, 1))


###################################################
### code chunk number 7: FME.Rnw:443-446
###################################################
DataLogV <- cbind(time = out$time,
            logV = out$logV + rnorm(sd = 0.45, n = length(out$logV)),
            sd = 0.45)


###################################################
### code chunk number 8: FME.Rnw:452-457
###################################################
ii    <- which (out$time %in% seq(0, 56, by = 4))
DataT <- cbind(time = out$time[ii],
               T = out$T[ii] + rnorm(sd = 4.5, n = length(ii)),
               sd = 4.5)
head(DataT)


###################################################
### code chunk number 9: FME.Rnw:511-516
###################################################
HIVcost <- function (pars) {
  out <- HIV(pars)
  cost <- modCost(model = out, obs = DataLogV, err = "sd")
  return(modCost(model = out, obs = DataT, err = "sd", cost = cost))
}


###################################################
### code chunk number 10: FME.Rnw:522-523
###################################################
HIVcost(pars)$model


###################################################
### code chunk number 11: cost
###################################################
plot(HIVcost(pars), xlab="time")


###################################################
### code chunk number 12: costfig
###################################################
plot(HIVcost(pars), xlab="time", legpos="topright")


###################################################
### code chunk number 13: sf
###################################################
# Local sensitivity analysis of bet on logV - code is hidden in text
ref  <- HIV(pars)
pert <- HIV (pars*c(1.1,1,1,1,1,1))
ss   <- (pert$logV-ref$logV)/pert$logV
plot(ref$time,ref$logV,type="l",lwd=2,xlab="time",ylab="logV",
  main="local sensitivity, parameter bet")
lines(pert$time,pert$logV,lwd=2,lty=2)
arrseq <- seq(9,36,6)#c(10,30,50,70,90)
arrows(ref$time[arrseq],ref$logV[arrseq],
  ref$time[arrseq],pert$logV[arrseq], length= 0.075)

legend("bottomright",c("bet=2e-5","bet=2.2e-5"),lty=c(1,2))
par(new=TRUE)
par(fig=c(0.5,0.99,0.5,0.95))
plot(ref$time,ss,type="l",lwd=2,
   xlab="",ylab="",axes=FALSE,frame.plot=TRUE)
points(ref$time[arrseq],ss[arrseq])
title("Sensitivity functions ",line=0.5,cex.main=1)
par(fig=c(0,1,0,1))


###################################################
### code chunk number 14: sf
###################################################
# Local sensitivity analysis of bet on logV - code is hidden in text
ref  <- HIV(pars)
pert <- HIV (pars*c(1.1,1,1,1,1,1))
ss   <- (pert$logV-ref$logV)/pert$logV
plot(ref$time,ref$logV,type="l",lwd=2,xlab="time",ylab="logV",
  main="local sensitivity, parameter bet")
lines(pert$time,pert$logV,lwd=2,lty=2)
arrseq <- seq(9,36,6)#c(10,30,50,70,90)
arrows(ref$time[arrseq],ref$logV[arrseq],
  ref$time[arrseq],pert$logV[arrseq], length= 0.075)

legend("bottomright",c("bet=2e-5","bet=2.2e-5"),lty=c(1,2))
par(new=TRUE)
par(fig=c(0.5,0.99,0.5,0.95))
plot(ref$time,ss,type="l",lwd=2,
   xlab="",ylab="",axes=FALSE,frame.plot=TRUE)
points(ref$time[arrseq],ss[arrseq])
title("Sensitivity functions ",line=0.5,cex.main=1)
par(fig=c(0,1,0,1))


###################################################
### code chunk number 15: FME.Rnw:616-618
###################################################
Sfun <- sensFun(HIVcost, pars)
summary(Sfun)


###################################################
### code chunk number 16: sfun
###################################################
plot(Sfun, which = c("logV", "T"), xlab="time", lwd = 2)


###################################################
### code chunk number 17: sfunfig
###################################################
plot(Sfun, which = c("logV", "T"), xlab="time", lwd = 2)


###################################################
### code chunk number 18: sfunp
###################################################
pairs(Sfun, which = c("logV", "T"), col = c("blue", "green"))


###################################################
### code chunk number 19: sfunpfig
###################################################
pairs(Sfun, which = c("logV", "T"), col = c("blue", "green"))


###################################################
### code chunk number 20: FME.Rnw:749-751
###################################################
ident <- collin(Sfun)
head(ident, n = 20)


###################################################
### code chunk number 21: coll
###################################################
plot(ident, log = "y")


###################################################
### code chunk number 22: collfig
###################################################
plot(ident, log = "y")


###################################################
### code chunk number 23: FME.Rnw:783-784
###################################################
collin(Sfun, parset = c("bet", "rho", "delt", "c", "lam", "n"))


###################################################
### code chunk number 24: FME.Rnw:791-792
###################################################
collin(Sfun, N = 5)


###################################################
### code chunk number 25: FME.Rnw:800-802
###################################################
collin(Sfun, parset = c("bet", "rho", "delt", "c", "lam"), which = "logV")
collin(Sfun, parset = c("bet", "rho", "delt", "c", "lam"), which = "T")


###################################################
### code chunk number 26: FME.Rnw:834-836
###################################################
HIVcost2 <- function(lpars)
  HIVcost(c(exp(lpars), n = 900))


###################################################
### code chunk number 27: FME.Rnw:844-848
###################################################
Pars <- pars[1:5] * 2
Fit <- modFit(f = HIVcost2, p = log(Pars))
exp(coef(Fit))
deviance(Fit)


###################################################
### code chunk number 28: FME.Rnw:854-856
###################################################
ini   <- HIV(pars = c(Pars, n = 900))
final <- HIV(pars = c(exp(coef(Fit)), n = 900))


###################################################
### code chunk number 29: fit
###################################################
par(mfrow = c(1,2))
plot(DataLogV, xlab = "time", ylab = "logV", ylim = c(7, 11))
lines(ini$time, ini$logV, lty = 2)
lines(final$time, final$logV)
legend("topright", c("data", "initial", "fitted"),
   lty = c(NA,2,1), pch = c(1, NA, NA))
plot(DataT, xlab = "time", ylab = "T")
lines(ini$time, ini$T, lty = 2)
lines(final$time, final$T)
par(mfrow = c(1, 1))


###################################################
### code chunk number 30: fitfig
###################################################
par(mfrow = c(1,2))
plot(DataLogV, xlab = "time", ylab = "logV", ylim = c(7, 11))
lines(ini$time, ini$logV, lty = 2)
lines(final$time, final$logV)
legend("topright", c("data", "initial", "fitted"),
   lty = c(NA,2,1), pch = c(1, NA, NA))
plot(DataT, xlab = "time", ylab = "T")
lines(ini$time, ini$T, lty = 2)
lines(final$time, final$T)
par(mfrow = c(1, 1))


###################################################
### code chunk number 31: FME.Rnw:994-996
###################################################
var0 <- Fit$var_ms_unweighted
cov0 <- summary(Fit)$cov.scaled * 2.4^2/5


###################################################
### code chunk number 32: FME.Rnw:997-999 (eval = FALSE)
###################################################
## MCMC <- modMCMC(f = HIVcost2, p = Fit$par, niter = 5000, jump = cov0,
##                 var0 = var0, wvar0 = 0.1, updatecov = 50)


###################################################
### code chunk number 33: FME.Rnw:1000-1001 (eval = FALSE)
###################################################
## save(MCMC, file="mcmc.Rdata")


###################################################
### code chunk number 34: FME.Rnw:1004-1005
###################################################
load("mcmc.Rdata")


###################################################
### code chunk number 35: FME.Rnw:1014-1015
###################################################
options(width = 50)


###################################################
### code chunk number 36: FME.Rnw:1018-1020
###################################################
MCMC$pars <- exp(MCMC$pars)
summary(MCMC)


###################################################
### code chunk number 37: mcmc
###################################################
plot(MCMC, Full = TRUE)


###################################################
### code chunk number 38: mcmcfig
###################################################
par(mar=c(4, 4, 3, 1) + .1)
plot(MCMC, Full = TRUE)


###################################################
### code chunk number 39: mcmc2
###################################################
pairs(MCMC, nsample = 1000)


###################################################
### code chunk number 40: mcmcfig2
###################################################
pairs(MCMC, nsample = 1000,cex.labels=1.4,cex=0.7)


###################################################
### code chunk number 41: FME.Rnw:1094-1095
###################################################
sR <- sensRange(func = HIV, parms = pars, parInput = MCMC$par)


###################################################
### code chunk number 42: sr
###################################################
plot(summary(sR), xlab = "time")


###################################################
### code chunk number 43: srfig
###################################################
plot(summary(sR), xlab = "time")


###################################################
### code chunk number 44: FME.Rnw:1145-1151
###################################################
parRange <- cbind(min = 0.75 * pars, max = 1.25 * pars)
crlfun <- function (pars)  return(meanVirus = mean(HIV(pars)$V))

CRL <- modCRL(fun = crlfun, parRange = parRange, num = 500)

cor(CRL)[7, ]


###################################################
### code chunk number 45: crl
###################################################
plot(CRL, ylab = "number of virions", trace = TRUE)


###################################################
### code chunk number 46: crlfig
###################################################
plot(CRL, ylab = "number of virions", cex=0.5, trace = TRUE)


