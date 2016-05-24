##==============================================================================
##            Applications of the ccl4 model
##==============================================================================

#par(mfrow=c(2, 2))
library(FME)

##------------------------------------------------------------------------------
##                           the Model
##------------------------------------------------------------------------------

# the model is a FORTRAN-implemented example from deSolve

##------------------------------------
## parameter values
##------------------------------------


Pm <- c(
  ## Physiological parameters
  BW = 0.182,   # Body weight (kg)
  QP = 4.0  ,   # Alveolar ventilation rate (hr^-1)
  QC = 4.0  ,   # Cardiac output (hr^-1)
  VFC = 0.08,   # Fraction fat tissue (kg/(kg/BW))
  VLC = 0.04,   # Fraction liver tissue (kg/(kg/BW))
  VMC = 0.74,   # Fraction of muscle tissue (kg/(kg/BW))
  QFC = 0.05,   # Fractional blood flow to fat ((hr^-1)/QC
  QLC = 0.15,   # Fractional blood flow to liver ((hr^-1)/QC)
  QMC = 0.32,   # Fractional blood flow to muscle ((hr^-1)/QC)
  ## Chemical specific parameters for chemical
  PLA = 16.17,  # Liver/air partition coefficient
  PFA = 281.48, # Fat/air partition coefficient
  PMA = 13.3,   # Muscle/air partition coefficient
  PTA = 16.17,  # Viscera/air partition coefficient
  PB = 5.487,   # Blood/air partition coefficient
  MW = 153.8,   # Molecular weight (g/mol)
  VMAX = 0.04321671, # Max. velocity of metabolism (mg/hr) -calibrated
  KM = 0.4027255,    # Michaelis-Menten constant (mg/l) -calibrated
  ## Parameters for simulated experiment
  CONC = 1000,  # Inhaled concentration
  KL = 0.02,    # Loss rate from empty chamber /hr
  RATS = 1.0,   # Number of rats enclosed in chamber
  VCHC = 3.8    # Volume of closed chamber (l)
)

##------------------------------------
## the observations
##------------------------------------

Obs <- ccl4data[ccl4data$initconc==1000,]
plot(ChamberConc ~ time, data=Obs, xlab="Time (hours)",
         xlim=range(c(0, ccl4data$time)),
         ylab="Chamber Concentration (ppm)",
         log="y", main = "ccl4model")

Pm["CONC"] <-1000


##------------------------------------
## The model run
##------------------------------------
ccl4run <- function(pars, times=unique(Obs$time)) {
  ## state variables
  y <- c(AI = 21,   # total mass , mg
         AAM = 0,
         AT = 0,
         AF = 0,
         AL = 0,
         CLT = 0, ### area under the conc.-time curve in the liver
         AM = 0   ### the amount metabolized (AM)
   )

   VCH <- Pm[["VCHC"]] - Pm[["RATS"]]*Pm[["BW"]]
   AI0 <- VCH * Pm[["CONC"]]*Pm[["MW"]]/24450
   y["AI"] <- AI0

  ## run the model:
  as.data.frame(ccl4model(times, y, pars))
}

out <- ccl4run(Pm)
lines(out$time, out$CP, lwd=2)

##------------------------------------------------------------------------------
##       Local sensitivity analysis : sensitivity functions
##------------------------------------------------------------------------------

## 1. Sensitivity functions
## All parameters are sensitivity parameters, all variables selected
Sens <- sensFun(func=ccl4run, parms=Pm, varscale=1)

##------------------------------------
## univariate sensitivity
##------------------------------------
(SF<-summary(Sens))
plot(SF)


## select the ones with highest sensitivity
pselect <- names(Pm)[which (SF$L2>1.5)]

##------------------------------------
## bivariate sensitivity
##------------------------------------

pairs(Sens[, pselect])
mtext(outer=TRUE, side=3, line=-1.5,
      "Sensitivity functions", cex=1.5)

##------------------------------------
## multivariate sensitivity
##------------------------------------

Coll<-collin(Sens[, pselect])
head(Coll)
tail(Coll)

plot(Coll, log="y")
abline(h=20, col="red")

## 'identifiable parameter combinations'
Coll[which( Coll$collinearity < 20), ]

##------------------------------------------------------------------------------
##       Global sensitivity analysis : Sensitivity ranges
##------------------------------------------------------------------------------

## 1. Define parameter range 0.8 - 1.2 their value (sensitive pars only)
parRange <- data.frame(min=Pm[pselect]*0.8, max=Pm[pselect]*1.2)
rownames(parRange) <- names(Pm[pselect])

parRange
## sensitivity range for sensitivity variable CP
Sr <- summary(sensRange(func=ccl4run, parms=Pm,
                sensvar="CP", parRange=parRange[1, ], num=100))
plot(Sr, xlab="time, hour",
     ylab="Chamber Concentration (ppm)", main="Sensitivity BW")

## and so on for other parameters

##------------------------------------------------------------------------------
## Fit the model to data
##------------------------------------------------------------------------------
## fitted parameters
fitPAR <- c("VMAX", "CONC", "KM")

# 1. Define the model residuals

Residuals <- function(P) {
  Pm[fitPAR]<-P
  out <- as.data.frame(ccl4run(Pm))
  return(out$CP-Obs$ChamberConc)
}

## Are the parameters identifiable based on the data?   (CONC is NOT)
collin(sensFun(f=Residuals, parms= c(VMAX=0.04, CONC=1000, KM=0.4), map=NULL))


(MrqFit<-modFit(p=c(VMAX=0.04, CONC=1000, KM=0.4),
             lower=c(0, 0, 0), f=Residuals))
summary(MrqFit)

# run model with the optimized value:
 Pm[fitPAR]<-MrqFit$par
 fitted <- ccl4run(Pm)

plot(ChamberConc ~ time, data=Obs, xlab="Time (hours)",
         xlim=range(c(0, ccl4data$time)),
         ylab="Chamber Concentration (ppm)",
         log="y", main = "ccl4model-fitted")
lines(fitted$time, fitted$CP, lwd=2)

##------------------------------------------------------------------------------
## MCMC application
##------------------------------------------------------------------------------

## estimate of parameter covariances (to update parameters) and the model variance
(sP <- summary(MrqFit))
s2prior<-sP$modVariance

## adaptive MCMC
MCMC <- modMCMC(f=Residuals, p=MrqFit$par, niter=1000,
  var0=s2prior, wvar0=0.1, updatecov=10, lower=c(0, 0, 0), upper=c(2, 2000, 2))

plot(MCMC, Full=TRUE)
pairs(MCMC)
cor(MCMC$pars)
cov(MCMC$pars)
sP$cov.scaled

sR<-sensRange(parms=Pm, parInput=MCMC$pars, func=ccl4run)
plot(SUMM<-summary(sR),  which ="CP")
points(ChamberConc ~ time, data=Obs)

## DRAM
MCMC2 <- modMCMC(f=Residuals, p=MrqFit$par, jump=0.5*MrqFit$par, niter=1000, ntrydr=3,
  var0=s2prior, wvar0=1, updatecov=10, lower=c(0, 0, 0), upper=c(2, 2000, 2))

plot(MCMC2, Full=TRUE)
pairs(MCMC2)

plot(summary(sensRange(parms=Pm, parInput=MCMC2$pars, func=ccl4run)), which ="CP")


