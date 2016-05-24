##
## JAK-STAT pathway: Numeric implementation of steady states
##
  
# Load required packages
  library(ggplot2)
  library(deSolve)
  library(cOde)
  library(dMod)
  library(trust)
  library(parallel)
  
## Prepare equations ----------------------------------------
  setwd("~/dMod/testing/JakStat")
# Generate equations from csv  
  eq <- generateEquations(read.csv("topology.csv"))
# Define observables  
  observables <- c(tSTAT = "STAT + pSTAT + 2*pSTATdimer + off_tSTAT",
                   tpSTAT = "s_STAT*(pSTAT + 2*pSTATdimer) + off_tpSTAT")
# Define forcings and fixed states
  forcings <- "pEpoR"
# Add observables to ODE  
  eq <- addObservable(observables, eq)
  model <- generateModel(eq, forcings = forcings, fcontrol = "einspline")
  
## Get data --------------------------------------------------
  
# Read in data.frame and plot
  mydata <- read.csv("pnas_data_original.csv")
  plotData(list(Epo = mydata)) + expand_limits(y=0) + geom_line()
# Extract forcings (pEpoR) from data  
  mydataReceptor <- subset(mydata, name=="pEpoR")
  myforc <- mydataReceptor[,c("name", "time", "value")]
  myforc <- rbind(data.frame(name = "pEpoR", time = c(-20, -1), value = 0.1), myforc)
  mydata <- subset(mydata, name!="pEpoR")  
  
## Parameter transformation ------------------------------------
  
# Collect all inner parameters  
  innerpars <- getSymbols(c(eq, names(eq), observables, names(observables)), exclude = forcings)
  names(innerpars) <- innerpars
    
# Define transformations
  trafo1 <- replaceSymbols(names(observables), observables, innerpars) # Observable initial values
  midpars <- getSymbols(trafo1)
  
  trafo3 <- replaceSymbols(midpars, paste0("exp(log", midpars, ")"), trafo1) # log-transform
  outerpars <- getSymbols(trafo3)
  
# Generate parameter transformation function
  p1 <- P(trafo1)
  p2 <- Pi(eq[1:9], parameters = "STAT")
  p3 <- P(trafo3)
  p <- p1 %o% p2 %o% p3
# Initialize outer parameters  
  pini <- rep(0, length(outerpars))
  names(pini) <- outerpars
  fixed <- c(pEpoR = myforc[1, 3])
  
## Model prediction ----------------------------------------
  
# Collect times from data and augment by additional time points
  timesD <- sort(unique(c(-20, -1, mydata$time)))
  times <- seq(min(timesD), max(timesD), len=250)
# Generate model prediction function
  x <- Xs(model$func, model$extended, forcings = myforc)
# Evaluate model at initial parameters and plot
  out <- x(times, p(pini, fixed))
  plotPrediction(list(states=out))
  
  
  
## Fit data ------------------------------------------------

# Define objective function
  prior <- rep(0, length(outerpars)); names(prior) <- outerpars
  myfn <- function(pp, fixed=NULL, deriv=TRUE) 
    wrss(res(mydata, x(timesD, p(pp, fixed=fixed), deriv = deriv))) + constraintL2(c(pp, fixed), prior, 5)
  
# Fit the data by a trust region algorithm
  myfit <- trust(myfn, pini + rnorm(length(pini), 0, .1), rinit=1, rmax=10, iterlim=500, fixed = fixed)
# Evaluate model at best-fitting parameter values and plot
  prediction <- x(times, p(myfit$argument, fixed))
  plotCombined(list(states=prediction, input=long2wide(mydataReceptor)), list(states=mydata, input=mydataReceptor))
  
