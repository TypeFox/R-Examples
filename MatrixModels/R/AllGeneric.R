##' Updates the mean vector mu given the linear predictor
##' gamma. Evaluate the residuals and the weighted sum of squared
##' residuals.
##'
##' Note that the offset is added to the linear predictor before
##' calculating mu.
##' The sqrtXwt matrix can be updated but the sqrtrwt should not be in
##' that the weighted sum of squared residuals should be calculated
##' relative to fixed weights.  Reweighting is done in a separate call.
##' @title Update the fitted mean response
##' @param respM a response module
##' @param gamma the value of the linear predictor before adding the offset
##' @param ...
##' @return updated respM
setGeneric("updateMu", function(respM, gamma, ...)
           standardGeneric("updateMu"))

##' Update the weights, sqrtrwt and sqrtXwt
##' @title Update the residual and X weights
##' @param respM a response module
##' @param ...
##' @return updated response module
setGeneric("updateWts", function(respM, ...)
           standardGeneric("updateWts"))

if (FALSE) {                            # don't need this generic in R
##' Set new values of the coefficients.  Can be called with a single
##' vector argument and with a pair of vectors, representing a base and
##' an increment, plus a step factor.
##' @title set new values of the coefficients
##' @param predM a predictor module
##' @param base coefficient base value
##' @param incr increment
##' @param step step factor, defaults to 0 in which case incr is ignored
##' @param ...
##' @return predM
setGeneric("setCoef", function(predM, base, incr, step = 0, ...) standardGeneric("setCoef"))
}

##' Update any internal structures associated with sqrtXwt and the
##' weighted residuals.  The "V" matrix is evaluated from X using the
##' sqrtXwt matrix and a Vtr vector is calculated.
##' @title Reweight Prediction Module Structure Internals
##' @param predM a predictor module
##' @param sqrtXwt the sqrtXwt matrix
##' @param wtres the vector of weighted residuals
##' @param ...
##' @return updated predM
setGeneric("reweightPred", function(predM, sqrtXwt, wtres, ...)
           standardGeneric("reweightPred"))

if (FALSE) {                       # don't need this generic in R
##' Return the gamma vector
##' @title
##' @param predM a predictor module
##' @param ...
##' @return X %*% coef
setGeneric("gammaInc", function(predM, ...)
           standardGeneric("gammaInc"))
}

##' Solve for the coefficients, usually in the form of
##' coef <- solve(predM@fac, predM@Vtr, system = "A")

##' The squared length of the intermediate solution is attached as an
##' attribute of the returned value.
##' @title solve for the coefficients or coefficient increment
##' @param predM
##' @param ...
##' @return coefficient vector or increment
setGeneric("solveCoef", function(predM, ...)
           standardGeneric("solveCoef"))


##------------ all these should wander to  stats4  eventually: -----------------

## Make resid() into a reasonable S4 generic (still dispatching for S3):
setMethod("resid", "ANY", function(object, ...) residuals(object, ...))

## ditto for fitted.values() & coefficients():
setMethod("fitted.values", "ANY", function(object, ...) fitted(object, ...))
setMethod("coefficients",  "ANY", function(object, ...) coef  (object, ...))

