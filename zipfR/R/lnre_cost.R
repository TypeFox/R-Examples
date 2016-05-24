##
## internal: various cost functions for parameter estimation of LNRE models
##

## common signature: lnre.cost(model, spc, m.max=15)

## CHI-SQUARED COST FUNCTION (lnre.cost.chisq)
## compute standard chi-squared statistic for observed data wrt. model
##  - observations are V_1, ..., V_m.max and V+ = V - (V_1 + ... + V_m.max)
##  - statistic assumes independence of these random variables
##    (!= multivariate chi-squared used for g.o.f.)
##  - means and standard deviations of variables are obtained from model
##  - s.d. is clamped to >= 3 to avoid scaling up small differences
##    for which normal approximation is unsuitable
lnre.cost.chisq <- function (model, spc, m.max=15)
{
  N <- N(spc)
  
  E.Vm <- EVm(model, 1:m.max, N)        # expected values according to LNRE model
  E.V <- EV(model, N)
  E.Vplus <- E.V - sum(E.Vm)

  V.Vm <- VVm(model, 1:m.max, N)        # variances according to LNRE model
  V.V <- VV(model, N)
  V.Vplus <- V.V - sum(V.Vm)            # since we assume independence of the V_k

  V.Vm <- pmax(V.Vm, 9)                 # clamp variance to >= 9 (i.e., s.d. >= 3)
  V.Vplus <- max(V.Vplus, 9)


  O.Vm <- Vm(spc, 1:m.max)              # observed values
  O.V <- V(spc)
  O.Vplus <- O.V - sum(O.Vm)

  X2 <- sum( (O.Vm - E.Vm)^2 / V.Vm ) + (O.Vplus - E.Vplus)^2 / V.Vplus
  X2
}

## LINEAR COST FUNCTION (lnre.cost.linear)
## sum of absolute values of differences between observed and expected values
##  - for vocabulary size V and spectrum elements V_1, ... , V_m.max
lnre.cost.linear <- function (model, spc, m.max=15)
{
  N <- N(spc)
  
  E.Vm <- EVm(model, 1:m.max, N)
  E.V <- EV(model, N)

  O.Vm <- Vm(spc, 1:m.max)
  O.V <- V(spc)

  abs(E.V - O.V) + sum( abs(E.Vm - O.Vm) )
}

## SMOOTH LINEAR COST FUNCTION (lnre.cost.smooth.linear)
## smoothed version of linear cost function which avoids kink of abs(x) at x=0
##  - sum of sqrt(1 + (O-E)^2) (min. cost = 1, approximates abs(O-E) for larger differences)
##  - for vocabulary size V and spectrum elements V_1, ... , V_m.max
lnre.cost.smooth.linear <- function (model, spc, m.max=15)
{
  N <- N(spc)
  
  E.Vm <- EVm(model, 1:m.max, N)
  E.V <- EV(model, N)

  O.Vm <- Vm(spc, 1:m.max)
  O.V <- V(spc)

  sqrt(1 + (E.V - O.V)^2) + sum( sqrt(1 + (E.Vm - O.Vm)^2) )
}

## MEAN SQUARED ERROR (MSE) COST FUNCTION (lnre.cost.mse)
## sum (or average) of squared differences between observed and expected values
##  - for vocabulary size V and spectrum elements V_1, ... , V_m.max
lnre.cost.mse <- function (model, spc, m.max=15)
{
  N <- N(spc)
  
  E.Vm <- EVm(model, 1:m.max, N)
  E.V <- EV(model, N)

  O.Vm <- Vm(spc, 1:m.max)
  O.V <- V(spc)

  ( (E.V - O.V)^2 + sum((E.Vm - O.Vm)^2) ) / (1 + m.max)
}


