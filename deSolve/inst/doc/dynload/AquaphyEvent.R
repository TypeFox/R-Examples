#---------------------------------------------------------------------------
# A phytoplankton model with uncoupled carbon and nitrogen assimilation
# as a function of light and Dissolved Inorganic Nitrogen (DIN) concentration
#
# The example demonstrates how to use forcing functions in compiled code
#
# before trying this code, the FORTRAN program has to be compiled
# this can be done in R:
# system("R CMD SHLIB AquaphyForcing.f")
# do make sure that this file is in the working directory...
# (if not, use setwd() )
#---------------------------------------------------------------------------

library(deSolve)

##==============================================================================
## Running the aquaphy model with light and dilution as forcing functions...
##==============================================================================

parameters <- c(maxPhotoSynt   = 0.125,      # mol C/mol C/hr
                rMortPHY       = 0.001,      # /hr
                alpha          = -0.125/150, # uEinst/m2/s/hr
                pExudation     = 0.0,        # -
                maxProteinSynt = 0.136,      # mol C/mol C/hr
                ksDIN          = 1.0,        # mmol N/m3
                minpLMW        = 0.05,       # mol C/mol C
                maxpLMW        = 0.15,       # mol C/mol C
                minQuotum      = 0.075,      # mol C/mol C
                maxStorage     = 0.23,       # /h
                respirationRate= 0.0001,     # /h
                pResp          = 0.4,        # -
                catabolismRate = 0.06,       # /h
                dilutionRate   = 0.01,       # /h
                rNCProtein     = 0.2,        # mol N/mol C
                inputDIN       = 10.0,       # mmol N/m3
                rChlN          = 1,          # g Chl/mol N
                parMean        = 250.,       # umol Phot/m2/s
                dayLength      = 24.         # hours - 24 hrs light
                )

## =======================
## The initial conditions
## =======================
times <- seq(10, 24*20, 1)

state <- c(DIN    =  6.0,   # mmol N/m3
          PROTEIN = 20.0,   # mmol C/m3
          RESERVE =  5.0,   # mmol C/m3
          LMW     =  1.0)   # mmol C/m3

## ==================
## The events
## ==================

tevent <- seq(0,24*20, by=24)
le <- length(tevent)

eventdat <- data.frame(var="DIN",time = tevent, value=6, method="replace")

## ==================
## Running the model
## ==================

out <- aquaphy(times, state, parameters, events=list(data=eventdat))

## ======================
## Plotting model output
## ======================

par(oma = c(0, 0, 3, 0))

plot(out, which=c("PAR","Chlorophyll","DIN","NCratio"), xlab = "time, hours", 
  ylab = c("uEinst/m2/s","ug/l","mmolN/m3","molN/molC"), type="l", lwd=2)

mtext(outer = TRUE, side = 3, "AQUAPHY", cex = 1.5)

## =====================
## Summary model output
## =====================
t(summary(out))
