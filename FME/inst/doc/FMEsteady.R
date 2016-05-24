### R code from vignette source 'FMEsteady.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library("FME")
options(prompt = "> ")
options(width=70)
set.seed(357)


###################################################
### code chunk number 2: FMEsteady.Rnw:111-113
###################################################
par(mfrow=c(2, 2))
require(FME)


###################################################
### code chunk number 3: FMEsteady.Rnw:117-121
###################################################
pars <- c(upO2 = 360,  # concentration at upper boundary, mmolO2/m3
          cons = 80,   # consumption rate, mmolO2/m3/day
          ks = 1,      # O2 half-saturation ct, mmolO2/m3
          D = 1)       # diffusion coefficient, cm2/d


###################################################
### code chunk number 4: FMEsteady.Rnw:124-128
###################################################
n  <- 100                       # nr grid points
dx <- 0.05   #cm
dX <- c(dx/2, rep(dx, n-1), dx/2)  # dispersion distances; half dx near boundaries
X  <- seq(dx/2, len = n, by = dx)  # distance from upper interface at middle of box


###################################################
### code chunk number 5: FMEsteady.Rnw:134-152
###################################################
O2fun <- function(pars)
{
  derivs<-function(t, O2, pars)
  {
  with (as.list(pars),{

    Flux <- -D* diff(c(upO2, O2, O2[n]))/dX
    dO2  <- -diff(Flux)/dx - cons*O2/(O2 + ks)

    return(list(dO2, UpFlux = Flux[1], LowFlux = Flux[n+1]))
  })
 }

 # Solve the steady-state conditions of the model
 ox <- steady.1D(y = runif(n), func = derivs, parms = pars,
                 nspec = 1, positive = TRUE)
 data.frame(X = X, O2 = ox$y)
}


###################################################
### code chunk number 6: FMEsteady.Rnw:155-156
###################################################
ox <- O2fun(pars)


###################################################
### code chunk number 7: FMEsteady.Rnw:159-160
###################################################



###################################################
### code chunk number 8: O2plot
###################################################
plot(ox$O2, ox$X, ylim = rev(range(X)), xlab = "mmol/m3",
     main = "Oxygen", ylab = "depth, cm", type = "l", lwd = 2)


###################################################
### code chunk number 9: O2plot
###################################################
plot(ox$O2, ox$X, ylim = rev(range(X)), xlab = "mmol/m3",
     main = "Oxygen", ylab = "depth, cm", type = "l", lwd = 2)


###################################################
### code chunk number 10: FMEsteady.Rnw:181-185
###################################################
print(system.time(
Sens2 <- sensRange(parms = pars, func = O2fun, dist = "norm",
           num = 100, parMean = c(cons = 80), parCovar = 100)
))


###################################################
### code chunk number 11: sens
###################################################
par(mfrow = c(1, 2))
plot(Sens2, xyswap = TRUE, xlab = "O2",
     ylab = "depth, cm", main = "Sensitivity runs")
plot(summary(Sens2), xyswap = TRUE, xlab = "O2",
     ylab = "depth, cm", main = "Sensitivity ranges")
par(mfrow = c(1, 1))


###################################################
### code chunk number 12: sens
###################################################
par(mfrow = c(1, 2))
plot(Sens2, xyswap = TRUE, xlab = "O2",
     ylab = "depth, cm", main = "Sensitivity runs")
plot(summary(Sens2), xyswap = TRUE, xlab = "O2",
     ylab = "depth, cm", main = "Sensitivity ranges")
par(mfrow = c(1, 1))


###################################################
### code chunk number 13: FMEsteady.Rnw:210-211
###################################################
O2sens <- sensFun(func=O2fun,parms=pars)


###################################################
### code chunk number 14: FMEsteady.Rnw:216-217
###################################################
summary(O2sens)


###################################################
### code chunk number 15: pairs
###################################################
pairs(O2sens)


###################################################
### code chunk number 16: pairs
###################################################
pairs(O2sens)


###################################################
### code chunk number 17: FMEsteady.Rnw:235-236
###################################################
cor(O2sens[,-(1:2)])


###################################################
### code chunk number 18: FMEsteady.Rnw:241-243
###################################################
Coll <- collin(O2sens)
Coll


###################################################
### code chunk number 19: coll
###################################################
plot(Coll, log = "y")


###################################################
### code chunk number 20: coll
###################################################
plot(Coll, log = "y")


###################################################
### code chunk number 21: FMEsteady.Rnw:264-269
###################################################
O2dat <- data.frame(x = seq(0.1, 3.5, by = 0.1),
    y = c(279,260,256,220,200,203,189,179,165,140,138,127,116,
          109,92,87,78,72,62,55,49,43,35,32,27,20,15,15,10,8,5,3,2,1,0))
O2depth <- cbind(name = "O2", O2dat)        # oxygen versus depth
O2flux  <- c(UpFlux = 170)                  # measured flux


###################################################
### code chunk number 22: FMEsteady.Rnw:273-292
###################################################
O2fun2 <- function(pars)
{
  derivs<-function(t, O2, pars)
  {
  with (as.list(pars),{

    Flux <- -D*diff(c(upO2, O2, O2[n]))/dX
    dO2  <- -diff(Flux)/dx - cons*O2/(O2 + ks)

    return(list(dO2,UpFlux = Flux[1], LowFlux = Flux[n+1]))
    })
  }

 ox <- steady.1D(y = runif(n), func = derivs, parms = pars, nspec = 1,
                   positive = TRUE, rtol = 1e-8, atol = 1e-10)

 list(data.frame(x = X, O2 = ox$y),
      UpFlux = ox$UpFlux)
}


###################################################
### code chunk number 23: FMEsteady.Rnw:298-314
###################################################
Objective <- function (P)
{
 Pars <- pars
 Pars[names(P)]<-P
 modO2 <- O2fun2(Pars)

 # Model cost: first the oxygen profile
 Cost  <- modCost(obs = O2depth, model = modO2[[1]],
                  x = "x", y = "y")

 # then the flux
 modFl <- c(UpFlux = modO2$UpFlux)
 Cost  <- modCost(obs = O2flux, model = modFl, x = NULL, cost = Cost)

 return(Cost)
}


###################################################
### code chunk number 24: FMEsteady.Rnw:318-323
###################################################
print(system.time(
sF<-sensFun(Objective, parms = pars)
))
summary(sF)
collin(sF)


###################################################
### code chunk number 25: FMEsteady.Rnw:330-336
###################################################
collin(sF, parset = c("upO2", "cons", "ks"))
print(system.time(
Fit <- modFit(p = c(upO2 = 360, cons = 80, ks = 1),
                  f = Objective, lower = c(0, 0, 0))
                  ))
(SFit<-summary(Fit))


###################################################
### code chunk number 26: res
###################################################
plot(Objective(Fit$par), xlab = "depth", ylab = "",
       main = "residual", legpos = "top")


###################################################
### code chunk number 27: res
###################################################
plot(Objective(Fit$par), xlab = "depth", ylab = "",
       main = "residual", legpos = "top")


###################################################
### code chunk number 28: FMEsteady.Rnw:355-358
###################################################
Pars <- pars
Pars[names(Fit$par)] <- Fit$par
modO2 <- O2fun(Pars)


###################################################
### code chunk number 29: fit
###################################################
plot(O2depth$y, O2depth$x, ylim = rev(range(O2depth$x)), pch = 18,
     main = "Oxygen-fitted", xlab = "mmol/m3", ylab = "depth, cm")
lines(modO2$O2, modO2$X)


###################################################
### code chunk number 30: fit
###################################################
plot(O2depth$y, O2depth$x, ylim = rev(range(O2depth$x)), pch = 18,
     main = "Oxygen-fitted", xlab = "mmol/m3", ylab = "depth, cm")
lines(modO2$O2, modO2$X)


###################################################
### code chunk number 31: FMEsteady.Rnw:380-382
###################################################
Covar   <- SFit$cov.scaled * 2.4^2/3
s2prior <- SFit$modVariance


###################################################
### code chunk number 32: FMEsteady.Rnw:385-391
###################################################
print(system.time(
MCMC <- modMCMC(f = Objective, p = Fit$par, jump = Covar,
     niter = 1000, ntrydr = 2, var0 = s2prior, wvar0 = 1,
     updatecov = 100, lower = c(NA, NA, 0))
))
MCMC$count


###################################################
### code chunk number 33: mcmcplot
###################################################
plot(MCMC,Full=TRUE)


###################################################
### code chunk number 34: mcmcplot
###################################################
plot(MCMC,Full=TRUE)


###################################################
### code chunk number 35: mcmchist
###################################################
hist(MCMC, Full = TRUE)


###################################################
### code chunk number 36: mcmchist
###################################################
hist(MCMC, Full = TRUE)


###################################################
### code chunk number 37: mcmcpairs
###################################################
pairs(MCMC, Full = TRUE)


###################################################
### code chunk number 38: mcmcpairs
###################################################
pairs(MCMC, Full = TRUE)


###################################################
### code chunk number 39: FMEsteady.Rnw:437-439
###################################################
summary(MCMC)
cor(MCMC$pars)


###################################################
### code chunk number 40: mcmcran2
###################################################
plot(summary(sensRange(parms = pars, parInput = MCMC$par, f = O2fun, num = 500)),
  xyswap = TRUE)
points(O2depth$y, O2depth$x)


###################################################
### code chunk number 41: mcmcran2
###################################################
plot(summary(sensRange(parms = pars, parInput = MCMC$par, f = O2fun, num = 500)),
  xyswap = TRUE)
points(O2depth$y, O2depth$x)


