## =============================================================================
## A model of BOD + O2 in a river
##
## checked 24-08-09
## =============================================================================

## Biochemical Oxygen Demand (BOD) and oxygen (O2) dynamics
## in a river - example from function steady.1D in rootSolve

library(FME)
par(mfrow=c(2, 2))

# The observed data

Data<- matrix(nc=2, byrow=TRUE, data =c(
   50,   9.5,  1050,   3.5,  2050,   5.3,  3050,  10.5,
 4050,  35.7,  5050, 112.1,  6050, 180.6,  7050, 236.1,
 8050, 268.2,  9050, 284.8))
colnames(Data) <- c("dist","O2")

##------------------------------------------------------------------------------
##                           the Model
##------------------------------------------------------------------------------

#==================#
# Model equations
#==================#
O2BOD <- function (pars) {
  derivs <- function(t, state, pars)  {
    with(as.list(pars), {
      BOD <- state[1:N]
      O2  <- state[(N+1):(2*N)]
      ## BOD dynamics
      FluxBOD <-  v*c(BOD_0, BOD)  # fluxes due to water flow (advection)
      FluxO2  <-  v*c(O2_0, O2)
      BODrate <- r*BOD*O2/(O2+ks) # 1-st order consumption, Monod in oxygen

      ##rate of change = flux gradient - consumption  + reaeration (O2)
      dBOD         <- -diff(FluxBOD)/dx  - BODrate
      dO2          <- -diff(FluxO2)/dx   - BODrate + p*(O2sat-O2)

      return(list(c(dBOD=dBOD, dO2=dO2), BODrate=BODrate))
    })
  }    # END derivs

  ## solving the model to steady-state
  ## initial guess
  state <- c(rep(200, N), rep(200, N))

  ## steady-state solution
  out <-steady.1D (y=state, func=derivs, parms=pars,
           nspec=2, pos=TRUE, names=c("BOD", "O2"))

  data.frame(dist=x, unlist(out$y), BODrate=out$BODrate)
} # end BODO2

#==================#
# Model parameters
#==================#
pars <- c(
v       = 1e2,       # velocity, m/day
r       = 1,         # /day, first-order decay of BOD
p       = 1,         # /day, air-sea exchange rate
ks      = 1,         # mmol O2/m3, half-saturation conc
O2sat   = 300,       # mmol/m3 saturated oxygen conc
O2_0    = 50 ,       # mmol/m3 riverine oxygen conc
BOD_0   = 1500)      # mmol/m3 riverine BOD concentration

#==================#
# Model grid
#==================#
dx      = 100       # grid size, meters
x       <- seq(dx/2, 10000, by=dx)  # m, distance from river
N       <- length(x)

#==================#
# Model solution
#==================#

out <- O2BOD(pars)

# initial oxygen concentration
O2_in <- out$O2

plot(x, O2_in, xlab= "Distance from river",
     ylab="Oxygen, mmol/m3", main="O2-BOD model", type="l")

##------------------------------------------------------------------------------
##       Global sensitivity analysis : Sensitivity ranges        ##
##------------------------------------------------------------------------------

# 1. Sensitivity parameter ranges
parRange=matrix(nc=2, data=c(500, 600, 0.1, 0.5, 0.5, 2), byrow=TRUE)
rownames(parRange) <- c("v", "r", "p")

# 2. Calculate sensitivity ranges for O2, BOD, and BODrate

print(system.time(
  Sens<-summary(SS<-sensRange(parms=pars, map=1,
                func=O2BOD, dist="unif", parRange=parRange, num=100))
))

# 3. Plot the results...
plot(Sens)

##------------------------------------------------------------------------------
##       Global sensitivity analysis : What-if scenarios
##------------------------------------------------------------------------------

## The effect of river flow on the oxygen concentration at the river outflow.

## 1. Define Sensitivity parameter range
parRange <- data.frame(min=10, max=5000)
rownames(parRange) <- "v"
parRange

crlmodel <- function (pars) {
  ## steady-state solution
  out   <- O2BOD(pars)
  c(O2end=out$O2[100], BODend=out$BOD[100], MeanBODrate = mean(out$BODrate))
}

## there is no mapping variable here...
crl <-modCRL(parms=pars, func=crlmodel, dist="grid",
                parRange=parRange, num=100)

head(crl)
plot(crl, xlab="velocity,m/day")

##------------------------------------------------------------------------------
##       Local sensitivity analysis : sensitivity functions
##------------------------------------------------------------------------------

## All parameters, all variables...
Sens <- sensFun(parms=pars, func=O2BOD)

## univariate sensitivity
summary(Sens)

## bivariate sensitivity
pairs(Sens)
pairs(Sens, which=c("O2", "BOD", "BODrate"))

mtext(outer=TRUE, side=3, line=-1.5,
      "Sensitivity functions", cex=1.5)

cor(Sens[, -(1:2)])

## multivariate sensitivity
Coll <- collin(Sens)
head(Coll)
tail(Coll)

plot(Coll, log="y")
abline(h=20, col="red")

##------------------------------------------------------------------------------
##                 Fitting the model to the data
##------------------------------------------------------------------------------

## three parameters are fitted: r, p, ks

## 1. Define a model cost function
Objective <- function(X) {
  pars[c("r", "p", "ks")]<-X                       # set parameter values
  out <- O2BOD(pars)
  modCost(model=out, obs=Data, x="dist")           # return SSR between model and data
}

## 2. Are these identifiable?  (collinearity of ~20 is critical)
collin(sensFun(func=Objective, parms=c(r=0.5, p=0.5, ks=1), varscale=1))

## 3. nlminb finds the minimum; parameters constrained to be > 0
print("Port")
print(system.time(Fit<-modFit(p=c(r=0.5, p=0.5, ks=1),
                  f=Objective, lower=c(0, 0, 0), method="Port")))
summary(Fit)

## 4. run the model with best-fit parameters
pars[c("r", "p", "ks")]<-Fit$par
O2new    <- O2BOD(pars)

## Plotting output
plot(Data, pch=18, cex=2, xlab= "Distance from river",
     ylab="mmol/m3", main="Oxygen", col="darkblue", ylim=c(0, 300))

lines(x, O2_in, col="darkgrey")
lines(O2new$dist, O2new$O2, lwd=2)
legend("bottomright", c("initial guess", "fitted"),
       col=c("darkgrey", "black"), lwd=c(1, 2))

##------------------------------------------------------------------------------
## Fit the model to data using the Levenberg-Marquardt algorithm
##------------------------------------------------------------------------------

print("Levenberg-Marquardt")
print(system.time(FitMrq <- modFit(p=c(r=0.5, p=.5, ks=1),
                 f=Objective, lower=c(0, 0, 0))))

summary(FitMrq)

Cost <- Objective(FitMrq$par)
plot(Cost, main="residuals")

##------------------------------------------------------------------------------
## MCMC
##------------------------------------------------------------------------------

## estimate of parameter covariances (to update parameters) and the model variance
sP <-summary(FitMrq)
Covar   <- sP$cov.scaled * 2.4^2/2
s2prior <- sP$modVariance

## adaptive Metropolis
MCMC <- modMCMC(f=Objective, p=FitMrq$par, jump=Covar, niter=1000,
                var0=s2prior, wvar0=0.2, updatecov=100)

plot(MCMC, Full=TRUE)
pairs(MCMC)
hist(MCMC, Full=TRUE)

cor(MCMC$pars)
cov(MCMC$pars)
sP$cov.scaled

sR<-sensRange(parms=pars, parInput=MCMC$pars, func=O2BOD)
plot(summary(sR))

