## Library dependencies and plot theme --------------------------------------------

library(deSolve)
library(bvpSolve)
library(trust)
library(parallel)
library(ggplot2)
library(ggthemes)
library(R2CdeSolve)

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

# Set list of inputs (part of the estimation)
inputs <- c(
  # "inp1", "inp2", "inp3", ...
  )

# Add observable ODEs to the original ODEs (choose one of the two options, or combine them)
# Alternatively, an observation function may be used, see Section "Prediction Function".
f <- variableTransformation(observables, f)
f <- addObservable(observables, f)

# Generate the model C files, compile them and return a list with func and extended.
model1 <- generateModelIE(f, observed = names(observables), inputs = inputs, forcings = forcings, compile = TRUE, nGridpoints = 101)

## Parameter Transformations -------------------------------------------

# Define inner parameters (parameters occurring in the equations except forcings)
innerpars <- getSymbols(c(f, names(f)))
if(!is.null(forcings)) innerpars <- innerpars[!(innerpars%in%forcings)]
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



## Data ----------------------------------------------------------------------

# The data sheet must contain data for each condition for all observables
# AND inputs

conditions <- c(
  #"condition1", "condition2", ...
  )

datasheet <- read.table("datafile.csv") # with columns condition, name, time, value, sigma

data <- lapply(conditions, function(mycondition) subset(datasheet, condition == mycondition))
names(data) <- conditions

## Prediction Function ------------------------------------------------------


# Set different forcings per condition
timesF <- seq(0, 100, by=0.1)
u1 <- data.frame(name = "u1", time = timesF, value = 1*dnorm(timesF, 0, 5))
u2 <- data.frame(name = "u2", time = timesF, value = 3*dnorm(timesF, 0, 5))
u3 <- data.frame(name = "u3", time = timesF, value = 8*dnorm(timesF, 0, 5))

# Set condition-specific parameter transformations and generate p2p function
trafo1 <- replaceSymbols("logp1", "logp1_1", trafo)
trafo2 <- replaceSymbols("logp1", "logp1_2", trafo)
trafo3 <- replaceSymbols("logp1", "logp1_3", trafo)

p1 <- P(trafo1)
p2 <- P(trafo2)
p3 <- P(trafo3)

# Generate prediction functions for the different conditions (depends on different forces 
# but not on different parameter transformations)
optionsData <- list(tau=0, ngrid=101, type="logspline")
x1 <- Xv(model1$func_fa, model1$func_l, u1, data[[1]], optionsBvp = list(), optionsOde = list(method="lsode"), optionsData = optionsData)
x2 <- Xv(model1$func_fa, model1$func_l, u2, data[[2]], optionsBvp = list(), optionsOde = list(method="lsode"), optionsData = optionsData)
x3 <- Xv(model1$func_fa, model1$func_l, u3, data[[3]], optionsBvp = list(), optionsOde = list(method="lsode"), optionsData = optionsData)

# Function for the total model prediction, returns a list of predictions 
x <- function(times, pouter, fixed=NULL, guess=NULL, ...) {
  
  l <- list(
    x1(times, p1(pouter, fixed), guess=guess[[1]], ...),
    x2(times, p2(pouter, fixed), guess=guess[[2]], ...),
    x3(times, p3(pouter, fixed), guess=guess[[3]], ...)
    # ...
    )
  names(l) <- conditions
  
  return(l)
  
}

## Objective functions ----------------------------

# Objective function for trustOptim()
# Attention: Do not use xlist as a variable name in the global environment
fn <- function(pOuter, fixed=NULL, ...) {
  
  n <- outerpars
  names(pOuter) <- outerpars
  
  xlist <<- x(times, pOuter, fixed, ...)
  cOuter <- constraint(pOuter, n, 0, 10)
  
  v <- Reduce("+", lapply(xlist, function(o) attr(o, "value"))) + cOuter$value
  
  return(as.numeric(v))
  
}

gr <- function(pOuter, fixed=NULL) {
  n <- outerpars
  names(pOuter) <- outerpars
  cOuter <- constraint(pOuter, n, 0, 10)
  g <- Reduce("+", lapply(xlist, function(o) attr(o, "gradient")[n])) + cOuter$gradient[n] 
  return(as.numeric(g))
}

