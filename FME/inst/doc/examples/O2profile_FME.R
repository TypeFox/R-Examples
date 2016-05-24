##==============================================================================
## A simple model of oxygen, diffusing along a spatial gradient,
## with imposed upper and lower boundary concentration
## oxygen is consumed at maximal fixed rate, monod limitation
##==============================================================================
par(mfrow=c(2, 2))
library(FME)


##------------------------------------------------------------------------------
##                           the Model
##------------------------------------------------------------------------------

# The model parameters
pars <- c(upO2=360,   # concentration at upper boundary, mmolO2/m3
          lowO2=10,   # concentration at lower boundary, mmolO2/m3
          cons=80,    # consumption rate, mmolO2/m3/day
          ks=1,       # O2 half-saturation ct, mmolO2/m3
          D=1)        # Diffusion coefficient
# The model grid
n  <- 100                       # nr grid points
dx <- 0.05   #cm
dX <- c(dx/2, rep(dx, n-1), dx/2)  # dispersion distances; half the grid size near boundaries
X  <- seq(dx/2, len=n, by=dx)     # distance from upper interface at middle of box


O2fun <- function(pars) {
  derivs<-function(t, O2, pars) {
    with (as.list(pars), {
      Flux <- -D* diff(c(upO2, O2, lowO2))/dX
      dO2  <- -diff(Flux)/dx-cons*O2/(O2+ks)
      return(list(dO2, UpFlux = Flux[1], LowFlux = Flux[n+1]))
    })
  }
  ## Solve the steady-state conditions of the model
  ox <- steady.band(y=runif(n), func=derivs, parms=pars, nspec=1, positive=TRUE)
  data.frame(X=X, O2=ox$y)
}

ox <- O2fun(pars)
# Plot results
plot(ox$O2, ox$X, ylim=rev(range(X)), xlab="mmol/m3",
     main="Oxygen", ylab="depth, cm", type="l", lwd=2)

##------------------------------------------------------------------------------
##       Global sensitivity analysis : Sensitivity ranges
##------------------------------------------------------------------------------

## 1. Sensitivity parameter range: consumption rate
pRange <- data.frame (min=60, max=100)
rownames(pRange) <- "cons"

## 2. Calculate sensitivity ranges for O2
##    model is solved 100 times, uniform parameter distribution (default)
Sens <- summary(sensRange(parms=pars, func=O2fun, num=100,
                parRange=pRange))

## same, now with normal distribution of consumption (mean = 80, variance=100)
Sens2 <- sensRange(parms=pars, func=O2fun, dist="norm",
           num=100, parMean=c(cons=80), parCovar=100)

## 3. Plot results
plot(summary(Sens2), xyswap=TRUE, xlab= "O2",
     ylab="depth, cm", main="Sensitivity O2 model")

##------------------------------------------------------------------------------
##       Local sensitivity analysis : sensitivity functions
##------------------------------------------------------------------------------

## Sensitivity functions
O2sens <- sensFun(func=O2fun, parms=pars)

## univariate sensitivity
summary(O2sens)
plot(O2sens)
plot(summary(O2sens))

## bivariate sensitivity
pairs(O2sens)

cor(O2sens[, -(1:2)])

## multivariate sensitivity
Coll <- collin(O2sens)
Coll
plot(Coll, log="y")

##------------------------------------------------------------------------------
##   Fitting the model to the data - using Levenberg-Marquardt
##------------------------------------------------------------------------------

O2fun2 <- function(pars) {
  derivs <- function(t, O2, pars) {
    with (as.list(pars), {
      Flux <- -D*diff(c(upO2, O2, lowO2))/dX
      dO2  <- -diff(Flux)/dx-cons*O2/(O2+ks)
      return(list(dO2, UpFlux = Flux[1], LowFlux = Flux[n+1]))
    })
  }
  ## Solve the steady-state conditions of the model
  ox <- steady.band(y=runif(n), func=derivs, parms=pars, nspec=1, positive=TRUE)
  ## return both the oxygen profile AND the fluxes at both ends
  list(data.frame(x=X, O2=ox$y), UpFlux=ox$UpFlux, LowFlux=ox$LowFlux)
}

## The data
O2dat <- data.frame(x=seq(0.1, 3.5, by=0.1),
    y = c(279, 260, 256, 220, 200, 203, 189, 179, 165, 140, 138, 127, 116,
          109, 92, 87, 78, 72, 62, 55, 49, 43, 35, 32, 27, 20, 15, 15, 10, 8, 5, 3, 2, 1, 0))
O2depth <- cbind(name="O2", O2dat)        # oxygen versus depth
O2flux  <- c(UpFlux=170, LowFlux=0)       # measured fluxes

## 1. Objective function to minimise; all parameters are fitted
Objective <- function (x) {
  Pars <- pars
  Pars[names(x)] <- x

  ## Solve the steady-state conditions of the model
  modO2 <- O2fun2(Pars)

  ## Model cost: first the oxygen profile
  Cost  <- modCost(obs=O2depth, model=modO2[[1]], x="x", y="y")

  ## then the fluxes
  modFl <- c(UpFlux=modO2$UpFlux, LowFlux=modO2$LowFlux)
  Cost  <- modCost(obs=O2flux, model=modFl, x=NULL, cost=Cost)

  return(Cost)
}

## 2. collinearity of the parameters
sF <- sensFun(Objective, parms=c(upO2=360, lowO2=10, cons=80, ks=1, D=1))
plot(sF)

summary(sF)
plot(summary(sF))

collin(sF)


## 3. find the minimum; parameters constrained to be > 0
print(system.time(Fit <- modFit(p=c(upO2=360, lowO2=10, cons=80, ks=1),
                  f=Objective, lower=c(0, 0, 0, 0))))
Fit
(SFit <- summary(Fit))

## 4. plot the residuals
plot(Objective(Fit$par), xlab="depth", ylab="", main="residual", legpos="top")

## 5. Show best-fit
BestPar <- pars
BestPar[names(Fit$par)] <- Fit$par

modO2 <- O2fun(BestPar)

plot(O2depth$y, O2depth$x, ylim=rev(range(O2depth$x)), pch=18,
     main="Oxygen-fitted", xlab="mmol/m3", ylab="depth, cm")
lines(modO2$O2, modO2$X)
Cost <- Objective(Fit$par)

##------------------------------------------------------------------------------
##  Run MCMC
##------------------------------------------------------------------------------

## 1. use parameter covariances of fit to update parameters
Covar   <- SFit$cov.scaled * 2.4^2/4

## mean squared residuals of fit = prior for model variance
s2prior <- Fit$ms

## adaptive Metropolis
MCMC <- modMCMC(f=Objective, p=Fit$par, jump=Covar, niter=1000, ntrydr=3,
                var0=s2prior, wvar0=1, updatecov=100, lower=c(NA, 0, NA, 0))

plot(MCMC, Full=TRUE)
hist(MCMC, Full=TRUE)

pairs(MCMC, Full=TRUE)
summary(MCMC)
cor(MCMC$pars)
plot(summary(sensRange(parms=pars, parInput=MCMC$par,
                       f=O2fun, num=100)), xyswap=TRUE)

## 2. mean variance of separate fitted variables are prior for model variance
## This does not work so well...
s2priorvar <- Fit$var_ms

## unless we artificially increase varaiance for low O2 flux
s2priorvar[2] <- 1

MCMC2 <- modMCMC(f=Objective, p=Fit$par, jump=Covar, niter=1000,
                var0=s2priorvar, wvar0=1, updatecov=10, lower=c(NA, 0, NA, 0))
plot(MCMC2, Full=TRUE)
hist(MCMC2, Full=TRUE)
pairs(MCMC2, Full=TRUE)
plot(summary(sensRange(parms=pars, parInput=MCMC2$par, f=O2fun, num=500)), xyswap=TRUE)

## 3. idem 2 but with delayed rejection
MCMC3 <- modMCMC(f=Objective, p=Fit$par, jump=Covar, niter=1000, ntrydr=3,
                var0=s2priorvar, wvar0=1, updatecov=10, lower=c(NA, 0, NA, 0))
plot(MCMC3, Full=TRUE)
hist(MCMC3, Full=TRUE)
pairs(MCMC3, Full=TRUE)

sR <- sensRange(parms=pars, func=O2fun, parInput=MCMC3$par, num=100)
plot(summary(sR), xyswap=TRUE)

points(O2depth$y, O2depth$x)
