
#'Log-Likelihood of a fitStMoMo object
#'
#'Returns the log-likelihood of a fitted Stochastic Mortality Model.
#'
#'@param object an object of class \code{fitStMoMo} representing a 
#'Stochastic Mortality Model fitted to some data.
#'@param ... other arguments.
#'
#'@return The log-likelihood of the fitted model.
#'@export
logLik.fitStMoMo <- function (object, ...) {
  structure(object$loglik, df = object$npar, nobs=object$nobs, class = "logLik")
}