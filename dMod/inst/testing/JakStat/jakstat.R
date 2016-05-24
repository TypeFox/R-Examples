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
  reactions <- read.csv("topology.csv")
  states <- colnames(reactions)[-(1:2)]
  volumes <- structure(rep(1, length(states)), names = states)
  volumes[-(1:3)] <- "Vnuc"
  eq <- generateEquations(reactions, volumes = volumes)
  #eq <- mergeReactions(c(1, 1, 2), eq)
  #eq <- mergeReactions(2:3, eq)
  
# Define observables  
  observables <- c(tSTAT = "s_tSTAT*(STAT + pSTAT + 2*pSTATdimer) + off_tSTAT",
                   tpSTAT = "s_tpSTAT*(pSTAT + 2*pSTATdimer) + off_tpSTAT")
# Define forcings and fixed states
  forcings <- "pEpoR"
  fixed.zero <- c("pSTAT", "pSTATdimer", "npSTATdimer", "nSTAT1", "nSTAT2", "nSTAT3", "nSTAT4", "nSTAT5")#"off_tpSTAT")
  fixed.one <- c("Vnuc", "STAT")
# Add observables to ODE  
  eq <- addObservable(observables, eq)
  model <- generateModel(eq, forcings = forcings, fixed = c(fixed.zero, fixed.one))
  
## Get data --------------------------------------------------
  
# Read in data.frame and plot
  mydata <- read.csv("pnas_data_original.csv")
  plotData(list(Epo = mydata)) + expand_limits(y=0) + geom_line()
# Extract forcings (pEpoR) from data  
  mydataReceptor <- subset(mydata, name=="pEpoR")
  myforc <- mydataReceptor[,c("name", "time", "value")]
  mydata <- subset(mydata, name!="pEpoR")  
  
## Parameter transformation ------------------------------------
  
# Collect all inner parameters  
  innerpars <- getSymbols(c(eq, names(eq), observables, names(observables)), exclude=forcings)
  names(innerpars) <- innerpars
  
# Define transformations
  trafo <- innerpars # Initialize parameter transformation by the identity
  trafo <- replaceSymbols(names(observables), observables, trafo) # Observable initial values
  trafo <- replaceSymbols(fixed.zero, "0", trafo) # Steady states / initial values
  trafo <- replaceSymbols(fixed.one, "1", trafo) # Steady states / initial values
  trafo <- replaceSymbols(innerpars, paste0("exp(log", innerpars, ")"), trafo) # log-transform
  
# Get outer parameters from trafo
  outerpars <- getSymbols(trafo)
  
# Generate parameter transformation function
  p <- P(trafo)
  plambda <- P(c(lambda = "exp(loglambda)"), parameters = c(outerpars, "loglambda"))
  
# Initialize outer parameters  
  pini <- rep(0, length(outerpars))
  names(pini) <- outerpars
  
## Model prediction ----------------------------------------
  
# Collect times from data and augment by additional time points
  timesD <- sort(unique(mydata$time))
  times <- seq(min(timesD), max(timesD), len=250)
# Generate model prediction function
  x <- Xs(model$func, model$extended, forcings = myforc)
# Evaluate model at initial parameters and plot
  out <- x(times, p(pini))
  plotPrediction(list(states=out))
  
  
  
## Fit data ------------------------------------------------

# Define objective function
  prior <- rep(0, length(outerpars)); names(prior) <- outerpars
  myfn <- function(pp, fixed=NULL, deriv=TRUE) 
    wrss(res(mydata, x(timesD, p(pp, fixed=fixed), deriv = deriv))) + priorL2(plambda(pp, fixed), prior, lambda = "lambda")
  
# Fit the data by a trust region algorithm
  fixed <- NULL
  pini <- pini[!names(pini)%in%names(fixed)]
  myfit <- trust(myfn, pini + rnorm(length(pini), 0, .1), rinit=1, rmax=10, iterlim=500, fixed = c(loglambda = 0, fixed))
# Evaluate model at best-fitting parameter values and plot
  prediction <- x(times, p(myfit$argument, fixed))
  print(plotCombined(list(states=prediction, input=long2wide(mydataReceptor)), list(states=mydata, input=mydataReceptor)))

# Compute profile for loglambda
  proflist.approx <- profile.trust(myfn, c(myfit$argument, loglambda = 0), "loglambda", limits=c(-20, 0), algoControl = list(gamma = 1, reg = .1), fixed = fixed)
  proflist.exact  <- profile.trust(myfn, c(myfit$argument, loglambda = 0), "loglambda", limits=c(-20, 0), algoControl = list(gamma = 0, W = "identity", reoptimize = TRUE), optControl = list(iterlim = 5), fixed = fixed)
  
  print(plotPaths(proflist.approx, proflist.exact))
  
  
  
  
# Compute profiles    
  myfit <- trust(myfn, pini, rinit = 1, rmax = 10, fixed = c(loglambda = -5, fixed), fterm = 0)
  bestfit <- myfit$argument
  
  proflist.approx <- do.call(c, mclapply(names(bestfit), function(n) profile.trust(myfn, bestfit, n, limits=c(-5, 5), algoControl = list(gamma = 1, reg = 1e-1), fixed = c(loglambda = -5, fixed)), mc.cores=length(bestfit)))
  proflist.exact  <- do.call(c, mclapply(names(bestfit), function(n) profile.trust(myfn, bestfit, n, limits=c(-5, 5), algoControl = list(gamma = 1e-2, reg = 1e-1, reoptimize = TRUE), optControl = list(iterlim = 5), fixed = c(loglambda = -5, fixed)), mc.cores=length(bestfit)))
  
  print(plotProfile(proflist.approx, proflist.exact))
  