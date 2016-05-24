## Extract Coefficients
## coef.glarma <- function(object, type = "both",  ...) {
##   switch(type, "beta" = object$delta[1:ncol(object$X)],
##          "ARMA" = object$delta[(ncol(object$X) + 1):length(object$delta)],
##          "both" =
##          list("ARMA" = object$delta[(ncol(object$X) + 1):length(object$delta)],
##               "beta" = object$delta[1:ncol(object$X)]))
## }

coef.glarma <- function(object, types = "all",  ...) {
  if (object$type == "NegBin"){
    switch(types, "beta" = object$delta[1:ncol(object$X)],
           "ARMA" = object$delta[(ncol(object$X) + 1):
                                 (ncol(object$X) +
                                  length(object$thetaLags) +
                                 length(object$phiLags))],
           "NB" = object$delta[length(object$delta)],
           "all" = list("ARMA" = object$delta[(ncol(object$X) + 1):
                            (ncol(object$X) +
                             length(object$thetaLags) +
                             length(object$phiLags))],
                        "NB" = object$delta[length(object$delta)],
                        "beta" = object$delta[1:ncol(object$X)]))

  }
  else{
    if(types == "NB") {
        stop(paste("Negative binomial parameter is not available in",
                   ifelse(object$type == "Poi", "Poisson distribution",
                                 "binomial distribution")))
    }

    switch(types, "beta" = object$delta[1:ncol(object$X)],
           "ARMA" = object$delta[(ncol(object$X) + 1):length(object$delta)],
           "all" =
             list("ARMA" = object$delta[(ncol(object$X) + 1):
                                         length(object$delta)],
                  "beta" = object$delta[1:ncol(object$X)]))
  }
}

## Extract Residuals
residuals.glarma <- function(object, ...) object$residuals
## resid.glarma <- function(object, ...) object$residuals

## Extract Fitted Values
fitted.glarma <- function(object,  ...) {
    if (!inherits(object, "glarma"))
        stop("method is only for glarma objects")
    return(object$fitted.values)
}

## fitted.glarma <- function(object, fit = "glarma", log = FALSE,  ...) {
##   if(fit == "glarma") {
##     if (log == FALSE) fits <- object$mu
##     else fits <- object$W
##   }
##   if(fit == "eta"){
##     if (log == FALSE) fits <- exp(object$eta)
##     else fits <- object$eta
##   }
##   fits
## }

## Extract AIC
extractAIC.glarma <- function(fit, ...) fit$aic

## Extract Log Likelihood
logLik.glarma <- function(object, deriv = 0, ...) {
  switch(as.character(deriv),
         "0" = object$logLik,
         "1" = object$logLikDeriv,
         "2" = object$logLikDeriv2)
}

## Extract Number of Observations
nobs.glarma <- function(object, ...) length(object$y)

## Extract Model Data Frame
model.frame.glarma <- function(formula, ...) {
    if (is.null(formula$offset)) {
        as.data.frame(cbind(y = formula$y, formula$X))
    } else {
        as.data.frame(cbind(y = formula$y, formula$X,
                            "(offset)" = formula$offset))
    }
}

