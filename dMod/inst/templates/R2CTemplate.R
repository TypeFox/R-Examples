## Library dependencies and plot theme --------------------------------------------

library(deSolve)
library(trust)
library(parallel)
library(ggplot2)
library(ggthemes)
library(cOde)
library(dMod)

ggplot <- function(...) ggplot2::ggplot(...) + theme_few() + scale_color_colorblind() + scale_fill_colorblind()
qplot <- function(...) ggplot2::qplot(...) + theme_few() + scale_color_colorblind() + scale_fill_colorblind()


## Model Definition ------------------------------------------------------

# Read in model csv
reactionlist <- read.csv("topology.csv") 

# Translate data.frame into equations
f <- generateEquations(reactionlist)

# Define new observables based on ODE states
observables <- c(
  # y1 = "s*x1 + off"
)

# Set list of forcings
forcings <- c(
  # "u1", "u2", "u3", ...
  )

# List of fixed parameters which are known beforehand
fixed <- c(
  # "fixed1", "fixed2", ...
  )

# Add observable ODEs to the original ODEs or use an observation function
# Choose one of the three options, or combine them
f <- variableTransformation(observables, f)
f <- addObservable(observables, f)
g <- Y(observables, f)


# Generate the model C files, compile them and return a list with func and extended.
model0 <- generateModel(f, fixed = fixed, forcings = forcings, jacobian = "inz.lsodes")

## Parameter Transformations -------------------------------------------

# Define inner parameters (parameters occurring in the equations except forcings)
# Add names(observables) if addObservables(observables, f) is used
innerpars <- getSymbols(c(f, names(f), observables), exclude=c(forcings, "time"))
names(innerpars) <- innerpars

# Define additional parameter constraints, e.g. steady-state conditions
# Parameters (left-hand side) are replaced in the right-hand side of consecutive lines by resolveRecurrence() 
constraints <- resolveRecurrence(c(
  #p1 = "p2 + p3",
  #p4 = "p1*p5"   
  ))

# Build up a parameter transformation (constraints, log-transform, etc.)
# Start with replacing initial value parameters of the observables
trafo <- replaceSymbols(names(observables), observables, innerpars)
# Then employ the other parameter constraints
trafo <- replaceSymbols(names(constraints), constraints, trafo)
# Then do a log-transform of all parameters (if defined as positive numbers)
trafo <- replaceSymbols(innerpars, paste0("exp(log", innerpars, ")"), trafo)


## Specify different conditions -----------------------------------------------------

conditions <- c(
  #"condition1", "condition2", ...
  )

# Set condition-specific parameter transformations and generate p2p function
trafoL <- lapply(conditions, function(con) trafo); names(trafoL) <- conditions

specific <- c("")
trafoL <- lapply(conditions, function(con) {
  replaceSymbols(specific, paste(specific, con, sep="_"), trafoL[[con]])
}); names(trafoL) <- conditions

pL <- lapply(conditions, function(con) P(trafoL[[con]])); names(pL) <- conditions


# Set different forcings per condition
timesF <- seq(0, 100, by=0.1)
uL <- list(
  data.frame(name = "u1", time = timesF, value = 1*dnorm(timesF, 0, 5)),
  data.frame(name = "u2", time = timesF, value = 3*dnorm(timesF, 0, 5)),
  data.frame(name = "u3", time = timesF, value = 8*dnorm(timesF, 0, 5))
); names(uL) <- conditions


# Specify prediction functions for the different conditions (depends on different forces 
# but not on different parameter transformations)

xL <- lapply(conditions, function(con) Xs(model0$func, model0$extended, uL[[con]])); names(xL) <- conditions

# Function for the total model prediction, returns a list of predictions (choose one of
# the two possibilities)
x <- function(times, pouter, fixed=NULL, ...) {
  
  out <- lapply(conditions, function(cond) xL[[cond]](times, pL[[cond]](pouter, fixed), ...))
  names(out) <- conditions
  return(out)
  
}

x <- function(times, pouter, fixed=NULL, ...) {
  
  out <- lapply(conditions, function(cond) {
    pinner <- pL[[cond]](pouter, fixed)
    prediction <- xL[[cond]](times, pinner, ...)
    observation <- g(prediction, pinner, attach = TRUE)
    return(observation)
  }); names(out) <- conditions
  return(out)
  
}


## Data ----------------------------------------------------------------------

datasheet <- read.table("datafile.csv") # with columns condition, name, time, value, sigma
data <- lapply(conditions, function(mycondition) subset(datasheet, condition == mycondition))
names(data) <- conditions

## Objective Functions -------------------------------------------------------

# Data times
timesD <- unique(sort(unlist(sapply(data, function(d) d$time))))

# Initalize parameters 
outerpars <- getSymbols(do.call(c, trafoL[conditions]))
prior <- rep(0, length(outerpars)); names(prior) <- outerpars
pouter <- rnorm(length(prior), prior, 1); names(pouter) <- outerpars

# Objective function for trust()
obj <- function(pouter, fixed=NULL, deriv=TRUE) {
  
  prediction <- x(timesD, pouter, fixed = fixed, deriv = deriv)
  out <- lapply(names(data), function(cn) wrss(res(data[[cn]], prediction[[cn]])))
  
  # Working with weak prior (helps avoiding runaway solutions of the optimization problem)
  cOuter <- constraintL2(pouter, prior, sigma = 10)
  
  cOuter + Reduce("+", out)
  
}


## Howto proceed -------------------------------------------------

# Predicting and plotting
times <- seq(min(timesD), max(timesD), len=100)
prediction <- x(times, pouter)
plotPrediction(prediction)
plotPrediction(prediction, name %in% names(observables))

# Fitting
plotData(data)
myfit <- trust(obj, pouter, rinit=1, rmax=10, iterlim=500)
prediction <- x(times, myfit$argument)
plotCombined(prediction, data)
plotCombined(prediction, data, name%in%names(observables))
plotCombined(prediction, data, name%in%names(observables)) + facet_grid(name~condition, scales="free")
plotObjective(myfit)

# Fitting from random positions
center <- pouter
sink("output.txt")
fitlist <- mclapply(1:100, function(i) {
  
  deviation <- rnorm(length(center), 0, 1)
  pars <- center + deviation
  
  out <- NULL
  myfit <- try(trust(obj, pars, rinit=1, rmax=10, iterlim=1000), silent=TRUE)
  if(!inherits(myfit, "try-error")) {
    out <- data.frame(index = i, 
                      chisquare = myfit$value, 
                      converged = myfit$converged, 
                      iterations = myfit$iterations, 
                      as.data.frame(as.list(myfit$argument)))
    cat("out", i, myfit$value, myfit$converged, myfit$iterations, "\n")
  }
  
  return(out)
  
}, mc.cores=24, mc.preschedule=FALSE)
sink()
fitlist <- do.call(rbind, fitlist[sapply(fitlist, class) == "data.frame"])
fitlist <- fitlist[order(fitlist$chisquare),]
save(fitlist, file="fitlist.rda")
bestfit <- unlist(fitlist[1,-(1:4)])
qplot(y = fitlist$chisquare)


# Profile likelihood
bestfit <- myfit$argument
profiles.approx <- do.call(c, mclapply(names(bestfit), function(n) profile.trust(obj, bestfit, n, limits=c(-3, 3), algoControl = list(gamma = 0)), mc.cores=4))
profiles.exact  <- do.call(c, mclapply(names(bestfit), function(n) profile.trust(obj, bestfit, n, limits=c(-3, 3), algoControl = list(gamma = 0, reoptimize = TRUE), optControl = list(iterlim = 10)), mc.cores=4))
plotProfile(profiles.approx, profiles.exact)
plotPaths(profiles.approx[1])
plotPaths(profiles.approx[c(1,3)])