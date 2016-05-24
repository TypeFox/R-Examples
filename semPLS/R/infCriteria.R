### Armin (first: 2012-09-21; Last: 2012-09-25)

### ic - information criteria
ic <- function(object, LV, criteria, ...){
  UseMethod("ic", object)
}

ic.sempls <- function(object, LV, criteria){
  criteriaOpt <- c("FPE", "AdjRsq", "AIC", "AICu", "AICc", "BIC",
                   "HQ", "HQc", "AIC2", "BIC2", "Cp", "GM")
  if(missing(criteria)) criteria <- criteriaOpt
  stopifnot(criteria %in% criteriaOpt)
  N <- object$N
  pred <- semPLS:::predecessors(object$model)[[LV]]
  pk <- length(pred)
  SSerrk <- deviance(object, LV)
  sapply(criteria, function(criteria)
         do.call(criteria, list(N, pk, SSerrk, object, LV)))
}

## Final Prediction Error
FPE <- function (N, pk, SSerrk, object, LV){
  SSerrk/N*(1 + (2*pk)/(N-pk))
}


## Adjusted R-sq; This value should match the one given by sempls package.
AdjRsq <- function (N, pk, SSerrk, object, LV){
  SStotal <- N-1
  1 - ((N-1)/(N-pk))*(SSerrk/SStotal)
}

## AIC
AIC <- function (N, pk, SSerrk, object, LV){
  N*(log(SSerrk/N) + (2*pk/N))
}

AIC2 <- function (N, pk, SSerrk, object, LV){
  ll <- 0.5 * (-N * (log(2 * pi) + 1 - log(N) + log(SSerrk)))
  -2 * ll + (2 * (pk + 1))
}


## Unbiased AIC
AICu <- function (N, pk, SSerrk, object, LV){
  N*(log(SSerrk/(N-pk)) + (2*pk/N))
}


## Corrected AIC
AICc <- function (N, pk, SSerrk, object, LV){
  N*(log(SSerrk/N) + (N+pk)/(N-pk-2))
}


## Bayesian IC
BIC <- function (N, pk, SSerrk, object, LV){
  N*(log(SSerrk/N) + pk*log(N)/N)
}

BIC2 <- function (N, pk, SSerrk, object, LV){
  ll <- 0.5 * (-N * (log(2 * pi) + 1 - log(N) + log(SSerrk)))
  -2 * ll + (pk + 1) * log(N)
}


## Hannan Quinn
HQ <- function (N, pk, SSerrk, object, LV){
  N*(log(SSerrk/N) + (2*pk*log(log(N)))/N)
}

## Corrected Hannan Quinn
HQc <- function (N, pk, SSerrk, object, LV){
  N*(log(SSerrk/N) + (2*pk*log(log(N)))/(N-pk-2))
}


### need saturated model:

## Mallows Cp
Cp <- function (N, pk, SSerrk, object, LV){
  objectFull <- fitsaturated(object, LV)
  MSerrFullModel <-  (1-rSquared(objectFull)[LV,]) # R-Squared from the saturated model.
  (SSerrk / MSerrFullModel) - (N - 2 * pk)
}

## Geweke Meese Criterion
GM <- function (N, pk, SSerrk, object, LV){
  objectFull <- fitsaturated(object, LV)
  MSerrFullModel <-  (1-rSquared(objectFull)[LV,]) # R-Squared from the saturated model. 
  (SSerrk/MSerrFullModel) + pk*log(N);
}
