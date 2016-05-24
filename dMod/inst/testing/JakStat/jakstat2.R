##
## JAK-STAT pathway
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
  reactions <- read.csv("topology_inputModel.csv")
  states <- colnames(reactions)[-(1:2)]
  volumes <- structure(rep(1, length(states)), names = states)
  volumes[-(1:5)] <- "Vnuc"
  eq <- generateEquations(reactions, volumes = volumes)
# Define observables
  observables <- c(tSTAT = "STAT + pSTAT + 2*pSTATdimer + offsettSTAT",
                   tpSTAT = "sSTAT*(pSTAT + 2*pSTATdimer + offsettpSTAT)",
                   pEpoR = "y1*y2")
  # Define forcings and fixed states
  forcings <- NULL
  zeroPars <- c("pSTAT", "pSTATdimer", "npSTATdimer", "nSTAT1", "nSTAT2", "nSTAT3", "nSTAT4", "nSTAT5", "y1")
  onePars <- c("y2")
# Add observables to ODE  
  eq <- addObservable(observables, eq)
  model <- generateModel(eq, fixed = c(zeroPars, onePars), nGridpoints = 101, 
                         jacobian = "inz.lsodes")
  
## Get data --------------------------------------------------
  
# Read in data.frame and plot
  mydata <- read.csv("pnas_data_original.csv")
  plotData(list(Epo = mydata)) + geom_line()

## Parameter transformation ------------------------------------
  
# Collect all inner parameters  
  innerpars <- getSymbols(c(eq, names(eq), observables, names(observables)), exclude=c(forcings, "time"))
  names(innerpars) <- innerpars
  
# Define transformations
  trafo <- replaceSymbols(names(observables), observables, innerpars) #observable initial values
  trafo <- replaceSymbols(zeroPars, "0", trafo) # steady states
  trafo <- replaceSymbols(onePars, "1", trafo) # steady states
  trafo <- replaceSymbols(innerpars, paste0("exp(", innerpars, ")"), trafo) #log-transform
  
# Get outer parameters from trafo
  outerpars <- getSymbols(trafo)
  
# Generate parameter transformation function
  p <- P(trafo)
  
# Initialize outer parameters  
  pini <- rep(0, length(outerpars))
  names(pini) <- outerpars
  
## Model prediction ----------------------------------------
  
# Collect times from data and augment by additional time points
  timesD <- sort(unique(mydata$time))
  times <- seq(min(timesD), max(timesD), len=250)
# Generate model prediction function
  x <- Xs(model$func, model$extended)
# Evaluate model at initial parameters and plot
  out <- x(times, p(pini))
  plotPrediction(list(states=out))
  
  
  
## Fit data ------------------------------------------------

# Define objective function
  prior <- rep(0, length(outerpars)); names(prior) <- outerpars
  myfn <- function(pp, fixed=NULL, deriv=TRUE) 
    wrss(res(mydata, x(timesD, p(pp, fixed=fixed), deriv = deriv))) + constraintL2(c(pp, fixed), prior, 20)
  
# Fit the data by a trust region algorithm
  myfit <- trust(myfn, pini + rnorm(length(pini), 0, .1), rinit=1, rmax=10, iterlim=500)
# Evaluate model at best-fitting parameter values and plot
  prediction <- x(times, p(myfit$argument))
  plotCombined(list(states=prediction, input=long2wide(mydataReceptor)), list(states=mydata, input=mydataReceptor))
  
## Multiple fits ------------------------------------------

# Fit repeatedly for different initial parameter values
# Return data.frame with the final objective value and 
# the best-fitting parameter vector
  fitlist <- mclapply(1:100, function(i) {
    
    mypini <- rnorm(length(outerpars), pini, 5)
    names(mypini) <- outerpars
    myfit <- trust(myfn, mypini, rinit=1, rmax=10, iterlim=100)
    out <- data.frame(chisquare = myfit$value, as.data.frame(as.list(myfit$argument)))
    
  }, mc.preschedule=FALSE, mc.cores=4)
# Sort the data.frame by increasing objective values and plot  
  fitlist <- do.call(rbind, fitlist[sapply(fitlist, class) == "data.frame"])
  fitlist <- fitlist[order(fitlist$chisquare),]
  qplot(y=fitlist[,1]) + scale_y_log10()

## Compute profiles of best fit -----------------------------
  
  
  bestfit <- unlist(fitlist[1,-1])
  proflist.approx <- do.call(c, mclapply(names(bestfit), function(n) profile.trust(myfn, bestfit, n, limits=c(-3, 3), algoControl = list(gamma = 1, reg = 0.1)), mc.cores=4))
  proflist.exact  <- do.call(c, mclapply(names(bestfit), function(n) profile.trust(myfn, bestfit, n, limits=c(-3, 3), algoControl = list(gamma = 0, reoptimize = TRUE), optControl = list(iterlim = 5)), mc.cores=4))
  
  plotProfile(proflist.approx, proflist.exact)
  