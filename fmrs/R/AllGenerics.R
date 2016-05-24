#' @title nobs method
#' @description Provides the number of observations in
#' an \code{\link{fmrsfit-class}}
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name nobs
#' @rdname nobs-methods
#' @param object An \code{\link{fmrsfit-class}}
#' @param ... Other possible arguments
#' @return An integer value
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' nObs = 500
#' dispersion = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrs.gendata(nObs = nObs, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), dispersion = dispersion,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrs.mle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDispersion = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' nobs(res.mle)
#' @exportMethod nobs
setGeneric("nobs", function(object, ...) standardGeneric("nobs"))
nobs.fmrsfit <- function(object, ...) {object@nobs}

#' @title ncov method
#' @description Provides the number of covariates of an FMRs model from
#' an \code{\link{fmrsfit-class}}
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name ncov
#' @rdname ncov-methods
#' @param object An \code{\link{fmrsfit-class}}
#' @param ... Other possible arguments
#' @return An integer value
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' nObs = 500
#' dispersion = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrs.gendata(nObs = nObs, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), dispersion = dispersion,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrs.mle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDispersion = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' ncov(res.mle)
#' @exportMethod ncov
setGeneric("ncov", function(object, ...) standardGeneric("ncov"))
ncov.fmrsfit <- function(object, ...) {object@ncov}

#' @title ncomp method
#' @description Provides the order of an FMRs model from
#' an \code{\link{fmrsfit-class}}
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name ncomp
#' @rdname ncomp-methods
#' @param object An \code{\link{fmrsfit-class}}
#' @param ... Other possible arguments
#' @return An integer value
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' nObs = 500
#' dispersion = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrs.gendata(nObs = nObs, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), dispersion = dispersion,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrs.mle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDispersion = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' ncomp(res.mle)
#' @exportMethod ncomp
setGeneric("ncomp", function(object, ...) standardGeneric("ncomp"))
ncomp.fmrsfit <- function(object, ...) {object@ncomp}

#' @title coefficients method
#' @description Provides the estimated regression coefficients from the
#' fitted FMRs model from an \code{\link{fmrsfit-class}}
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name coefficients
#' @rdname coefficients-methods
#' @param object An \code{\link{fmrsfit-class}}
#' @param ... Other possible arguments
#' @return A numeric array of dimension-\code{(nCov+1)}-\code{nComp}
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' nObs = 500
#' dispersion = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrs.gendata(nObs = nObs, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), dispersion = dispersion,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrs.mle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDispersion = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' coefficients(res.mle)
#' @exportMethod coefficients
setGeneric("coefficients",
           function(object, ...) standardGeneric("coefficients"))
coefficients.fmrsfit <- function(object, ...) {object@coefficients}

#' @title dispersion method
#' @description Provides the estimated dispersions of the fitted FMRs model
#' from from an \code{\link{fmrsfit-class}}
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name dispersion
#' @rdname dispersion-methods
#' @param object An \code{\link{fmrsfit-class}}
#' @param ... Other possible arguments
#' @return A numeric array of dimension-\code{(nCov+1)}-\code{nComp}
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' nObs = 500
#' dispersion = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrs.gendata(nObs = nObs, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), dispersion = dispersion,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrs.mle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDispersion = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' dispersion(res.mle)
#' @exportMethod dispersion
setGeneric("dispersion", function(object, ...) standardGeneric("dispersion"))
dispersion.fmrsfit <- function(object, ...) {object@dispersion}

#' @title mixProp method
#' @description Provides the estimated mixing proportions of an FMRs model
#'     form an \code{\link{fmrsfit-class}}
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name mixProp
#' @rdname mixProp-methods
#' @param object An \code{\link{fmrsfit-class}}
#' @param ... Other possible arguments
#' @return A numeric vector of length-\code{nComp}
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' nObs = 500
#' dispersion = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrs.gendata(nObs = nObs, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), dispersion = dispersion,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrs.mle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDispersion = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' mixProp(res.mle)
#' @exportMethod mixProp
setGeneric("mixProp", function(object, ...) standardGeneric("mixProp"))
mixProp.fmrsfit <- function(object, ...) {object@mixProp}

#' @title fitted method
#' @description Provides the fitted response of the fitted FMRs model from
#' an \code{\link{fmrsfit-class}}
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name fitted
#' @rdname fitted-methods
#' @param object An \code{\link{fmrsfit-class}}
#' @param ... Other possible arguments
#' @return A numeric array of dimension-\code{nObs}-\code{nComp}
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' nObs = 500
#' dispersion = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrs.gendata(nObs = nObs, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), dispersion = dispersion,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrs.mle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDispersion = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' head(fitted(res.mle))
#' @exportMethod fitted
setGeneric("fitted", function(object, ...) standardGeneric("fitted"))
fitted.fmrsfit <- function(object, ...) {object@fitted}

#' @title residuals method
#' @description Provides the residuals of the fitted FMRs model from
#' an \code{\link{fmrsfit-class}}
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name residuals
#' @rdname residuals-methods
#' @param object An \code{\link{fmrsfit-class}}
#' @param ... Other possible arguments
#' @return A numeric array of dimension-\code{nObs}-\code{nComp}
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' nObs = 500
#' dispersion = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrs.gendata(nObs = nObs, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), dispersion = dispersion,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrs.mle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDispersion = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' head(residuals(res.mle))
#' @exportMethod residuals
setGeneric("residuals", function(object, ...) standardGeneric("residuals"))
residuals.fmrsfit <- function(object, ...) {object@residuals}

#' @title weights method
#' @description Provides the weights of fitted observations for
#' each observation under all components of an FMRs model
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name weights
#' @rdname weights-methods
#' @param object An \code{\link{fmrsfit-class}}
#' @param ... Other possible arguments
#' @return A numeric array of dimension-\code{nObs}-\code{nComp}
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' nObs = 500
#' dispersion = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrs.gendata(nObs = nObs, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), dispersion = dispersion,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrs.mle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDispersion = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' head(weights(res.mle))
#' @exportMethod weights
setGeneric("weights", function(object, ...) standardGeneric("weights"))
weights.fmrsfit <- function(object, ...) {object@weights}

#' @title logLik method
#' @description Provides the estimated logLikelihood of an FMRs model from
#'     an \code{\link{fmrsfit-class}}
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name logLik
#' @rdname logLik-methods
#' @param object An \code{\link{fmrsfit-class}}
#' @param ... Other possible arguments
#' @return A numeric value
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' nObs = 500
#' dispersion = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrs.gendata(nObs = nObs, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), dispersion = dispersion,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrs.mle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDispersion = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' logLik(res.mle)
#' @exportMethod logLik
setGeneric("logLik", function(object, ...) standardGeneric("logLik"))
logLik.fmrsfit <- function(object, ...) {object@logLik}

#' @title BIC method
#' @description Provides the estimated BIC of an FMRs model from
#' an \code{\link{fmrsfit-class}}
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name BIC
#' @rdname BIC-methods
#' @param object An \code{\link{fmrsfit-class}}
#' @param ... Other possible arguments
#' @return A numeric value
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' nObs = 500
#' dispersion = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrs.gendata(nObs = nObs, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), dispersion = dispersion,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrs.mle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDispersion = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' BIC(res.mle)
#' @exportMethod BIC
setGeneric("BIC", function(object, ...) standardGeneric("BIC"))
BIC.fmrsfit <- function(object, ...) {object@BIC}

#' @title summary method
#' @description Displays the fitted FMRs model by showing the estimated
#' coefficients, dispersions and mixing proportions
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name summary
#' @rdname summary-methods
#' @param object An \code{\link{fmrsfit-class}}
#' @param ... Other possible arguments
#' @return Summary of the fitted FMRs model
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' nObs = 500
#' dispersion = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrs.gendata(nObs = nObs, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), dispersion = dispersion,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrs.mle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDispersion = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' summary(res.mle)
#' @exportMethod summary
setGeneric("summary", function(object, ...) standardGeneric("summary"))
summary.fmrsfit <- function(object, ...) {
  if(object@model == "FMR") {
    modelfmr = "Finite Mixture of Regression Models"
  }else if(object@disFamily == "lnorm"){
    modelfmr = "Finite Mixture of Accelerated Failure Time Regression
    Models \n  Log-Normal Sub-Distributions"
  }else{
    modelfmr = "Finite Mixture of Accelerated Failure Time Regression
    Models \n  Weibull Sub-Distributions"
  }
  cat("-------------------------------------------\n")
  cat("Fitted Model: \n")
  cat("-------------------------------------------\n")
  cat(" ", modelfmr, "\n")
  cat("  ", object@ncomp, " Components; ",
      object@ncov," Covariates; ", object@nobs,
      " samples.", sep = "")
  cat("\n\nCoefficients:\n")
  print.default(object@coefficients)
  cat("\nDispersions:\n")
  print.default(object@dispersion)
  cat("\nMixing Proportions:\n")
  print.default(object@mixProp)
  cat("\nLogLik: ", object@logLik, "; BIC: ", object@BIC, sep="")
  }

#' @title show method
#' @description Provides information about the fitted FMRs model
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @name show
#' @rdname show-methods
#' @param object An \code{\link{fmrsfit-class}}
#' @return Information about the fitted FMRs model
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' nObs = 500
#' dispersion = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrs.gendata(nObs = nObs, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), dispersion = dispersion,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrs.mle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDispersion = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#' show(res.mle)
#' @exportMethod show
setGeneric("show")

show.fmrsfit <- function(object) {
  if(object@model == "FMR") {
    modelfmr = "Finite Mixture of Regression Models"
  }else if(object@disFamily == "lnorm"){
    modelfmr = "Finite Mixture of Accelerated Failure Time Regression Models
    Log-Normal Sub-Distributions"
  }else{
    modelfmr = "Finite Mixture of Accelerated Failure Time Regression Models
    Weibull Sub-Distributions"
  }
  cat("An object of class '", class(object), "'\n", sep = "")
  cat(" ", modelfmr, "\n")
  cat("  ", object@ncomp, " Components; ",
      object@ncov," Covariates; ", object@nobs,
      " samples.\n", sep = "")
  }


#' @title  fmrs.mle method
#' @description Provides MLE for Finite Mixture of
#'     Accelerated Failure Time Regression Models or Finite Mixture of
#'     Regression Models. It also provides Ridge Regression.
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @family lnorm, norm, weibull
#' @name fmrs.mle
#' @rdname fmrs.mle-methods
#' @param y Responses (observations)
#' @param x Design matrix (covariates)
#' @param delta Censoring indicator vector
#' @param nComp Order (Number of components) of mixture model
#' @param disFamily A sub-distribution family. The options
#'     are \code{"norm"} for FMR models, \code{"lnorm"} for mixture of AFT
#'     regression models with Log-Normal sub-distributions,\code{"weibull"}
#'     for mixture of AFT regression models with Weibull sub-distributions,
#' @param initCoeff Vector of initial values for regression coefficients
#' including intercepts
#' @param initDispersion Vector of initial values for standard deviations
#' @param initmixProp Vector of initial values for proportion of components
#' @param lambRidge A positive value for tuning parameter in Ridge
#'     Regression or Elastic Net
#' @param nIterEM Maximum number of iterations for EM algorithm
#' @param nIterNR Maximum number of iterations for Newton-Raphson algorithm
#' @param conveps A positive value for avoiding NaN in computing divisions
#' @param convepsEM A positive value for treshold of convergence in
#'     EM algorithm
#' @param convepsNR A positive value for treshold of convergence in
#'     NR algorithm
#' @param porNR A positive interger for maximum number of searches in
#' NR algorithm
#' @param ... Other possible options
#' @keywords FMRs AFT Censored EM NR Ridge
#' @concept fmr, aft, mle, ridge, fmrs
#' @details Finite mixture of AFT regression models are represented as
#'     follows. Let \eqn{X} be the survival time with non-negative values,
#'     and \eqn{\boldsymbol{z} =(z_{1}, \ldots, z_{d})^{\top}}
#'     be a \eqn{d}-dimensional vector of covariates that may have an effect
#'     on \eqn{X}.
#'     If the survival time is subject to right censoring, then the observed
#'     response time is \eqn{T=\min \{Y, C\}},
#'     where \eqn{Y=\log X}, \eqn{C} is  logarithm of  the censoring time
#'     and \eqn{\delta=I_{\{y<c\}}} is the censoring indicator.
#'     We say that \eqn{V=(T,\delta,\boldsymbol z)} follows a finite mixture
#'     of AFT regression models of order \eqn{K} if the
#'     conditional density of \eqn{(T,\delta)} given \eqn{\boldsymbol z}
#'     has the form \deqn{f(t,\delta;\boldsymbol{z},\boldsymbol\Psi)
#'     =\sum\limits_{k=1}^{K}\pi_{k}[f_Y(t;\theta_{k}(\boldsymbol z),
#'     \sigma_{k})]^{\delta}[S_Y(t;\theta_{k}(\boldsymbol z)
#'     ,\sigma_{k})]^{1-\delta}[f_{C}(t)]^{1-\delta}[S_{C}(t)]^{\delta}}
#'     where \eqn{f_Y(.)} and \eqn{S_Y(.)} are respectively the density and
#'     survival functions of \eqn{Y}, \eqn{f_C(.)} and \eqn{S_C(.)} are
#'     respectively the density and survival functions of \eqn{C}; and
#'     \eqn{{\theta}_{k}(\boldsymbol{z})=h(\beta_{0k}+\boldsymbol{z}^{\top}
#'     \boldsymbol\beta_{k})} for a known link function \eqn{h(.)},
#'     \eqn{\boldsymbol\Psi=(\pi_{1},\ldots,\pi_{K},\beta_{01},\ldots,
#'     \beta_{0K},\boldsymbol\beta_{1},
#'     \ldots,\boldsymbol\beta_{K},\sigma_{1},
#'     \ldots,\sigma_{K})^{\top}} with \eqn{\boldsymbol\beta_{k}=
#'     (\beta_{k1},\beta_{k2},\ldots,\beta_{kd})^{\top}}
#'     and \eqn{0<\pi_{k}<1}
#'     with \eqn{\sum_{k=1}^{K}\pi_{k}=1}.
#'     The log-likelihood of a sample of size $n$ is formed
#'     as \deqn{\ell_{n}(\boldsymbol\Psi) =
#'     \sum\limits_{i=1}^{n}\log\sum\limits_{k=1}^{K}\pi_{k}\left[f_Y(t_{i},
#'     \theta_{k}({\boldsymbol z}_{i}),\sigma_{k})  \right]^{\delta_{i}}
#'     \left[S_Y(t_{i},\theta_{k}({\boldsymbol z}_{i}),
#'     \sigma_{k})\right]^{1-\delta_{i}}.}
#'     Note that we assume the censoring distribution is non-informative and
#'     hence won't play any role in the estimation process. We use EM and
#'     Newton-Raphson algorithms in our method to find the maximizer of
#'     above Log-Likelihood.
#' @references Shokoohi, F., Khalili, A., Asgharian, M. and Lin, S.
#'     (2016 submitted) Variable Selection in Mixture of Survival Models
#' @return An \code{\link{fmrsfit-class}} that includes parameter
#'     estimates of the specified FMRs model
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' nObs = 500
#' dispersion = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrs.gendata(nObs = nObs, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), dispersion = dispersion,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrs.mle(y = dat$y, x = dat$x, delta = dat$delta,
#'                     nComp = nComp, disFamily = "lnorm",
#'                     initCoeff = rnorm(nComp*nCov+nComp),
#'                     initDispersion = rep(1, nComp),
#'                     initmixProp = rep(1/nComp, nComp))
#' summary(res.mle)
#' @exportMethod fmrs.mle
setGeneric("fmrs.mle",
           function(y,delta,x,nComp,...) standardGeneric("fmrs.mle"))

#' @title  fmrs.tunsel method
#' @description Provides component-wise tuning parameters using BIC for
#'     Finite Mixture of Accelerated Failure Time Regression Models
#'     and Finite Mixture of Regression Models.
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @family lnorm, norm, weibull
#' @name fmrs.tunsel
#' @rdname fmrs.tunsel-methods
#' @param y Responses (observations)
#' @param x Design matrix (covariates)
#' @param delta Censoring indicator vector
#' @param nComp Order (Number of components) of mixture model
#' @param disFamily A sub-distribution family. The options
#'     are \code{"norm"} for FMR models,
#'     \code{"lnorm"} for mixture of AFT regression models with Log-Normal
#'     sub-distributions, \code{"weibull"} for mixture of AFT regression
#'     models with Weibull sub-distributions,
#' @param initCoeff Vector of initial values for regression coefficients
#'     including intercepts
#' @param initDispersion Vector of initial values for standard deviations
#' @param initmixProp Vector of initial values for proportion of components
#' @param penFamily Penalty name that is used in variable selection method.
#'     The available options are  \code{"lasso"}, \code{"adplasso"},
#'     \code{"mcp"}, \code{"scad"}, \code{"sica"} and \code{"hard"}.
#' @param lambRidge A positive value for tuniing parameter in Ridge
#'     Regression or Elastic Net
#' @param nIterEM Maximum number of iterations for EM algorithm
#' @param nIterNR Maximum number of iterations for Newton-Raphson algorithm
#' @param conveps A positive value for avoiding NaN in computing divisions
#' @param convepsEM A positive value for treshold of convergence in
#'     EM algorithm
#' @param convepsNR A positive value for treshold of convergence in
#'     NR algorithm
#' @param porNR A positive interger for maximum number of searches in
#' NR algorithm
#' @param gamMixPor Proportion of mixing parameters in the penalty. The
#'     value must be in the interval [0,1]. If \code{gamMixPor = 0}, the
#'     penalty structure is no longer mixture.
#' @param ... Other possible options
#' @keywords FMRs AFT Censored Tuning Ridge Regression LASSO Adaptive MCP
#' SCAD SICA
#' @concept fmr, aft, lasso, adplasso, mcp, scad, sica, ridge
#' @details The maximizer of penalized Log-Likelihood depends on selecting a
#'     set of good tuning parameters which is a rather thorny issue. We
#'     choose a value in an equally spaced set of values in
#'     \eqn{(0, \lambda_{max})} for a pre-specified \eqn{\lambda_{max}}
#'     that maximize the component-wise
#'     BIC, \deqn{\hat\lambda_{k} ={argmax}_{\lambda_{k}}BIC_k(\lambda_{k})=
#'     {argmax}_{\lambda_{k}}\left\{\ell^{c}_{k, n}
#'     (\hat{\boldsymbol\Psi}_{\lambda_{k}, k}) -
#'     |d_{\lambda_{k},k}| \log (n)\right\},}
#'     where \eqn{d_{\lambda_{k},k}=\{j:\hat{\beta}_{\lambda_{k},kj}\neq 0,
#'     j=1,\ldots,d\}} is the active set  excluding the intercept
#'     and \eqn{|d_{\lambda_{k},k}|}
#'     is its size. This approach is much faster than using an \code{nComp}
#'     by \code{nComp} grid to select the set \eqn{\boldsymbol\lambda} to
#'     maximize the penallized Log-Likelihood.
#' @references Shokoohi, F., Khalili, A., Asgharian, M. and Lin, S.
#'     (2016 submitted) Variable Selection in Mixture of Survival Models
#' @return An \code{\link{fmrstunpar-class}} that includes
#'     component-wise tuning parameter estimates that can be used in
#'     variable selection procedure.
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' nObs = 500
#' dispersion = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrs.gendata(nObs = nObs, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), dispersion = dispersion,
#'                      mixProp = mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrs.mle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDispersion = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#'
#' res.lam <- fmrs.tunsel(y = dat$y, x = dat$x, delta = dat$delta,
#'                       nComp = nComp, disFamily = "lnorm",
#'                       initCoeff = c(coefficients(res.mle)),
#'                       initDispersion = dispersion(res.mle),
#'                       initmixProp = mixProp(res.mle),
#'                       penFamily = "adplasso")
#' show(res.lam)
#' @exportMethod fmrs.tunsel
setGeneric("fmrs.tunsel",
           function(y,delta,x,nComp,...) standardGeneric("fmrs.tunsel"))

#' @title fmrs.varsel method
#' @description Provides variable selection and penalized MLE for
#'     Finite Mixture of Accelerated Failure Time Regression (FMAFTR) Models
#'     and Finite Mixture of Regression (FMR) Models.
#'     It also provide Ridge Regression and Elastic Net.
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @family lnorm, norm, weibull
#' @name fmrs.varsel
#' @rdname fmrs.varsel-methods
#' @param y Responses (observations)
#' @param x Design matrix (covariates)
#' @param delta Censoring indicators
#' @param nComp Order (Number of components) of mixture model
#' @param disFamily A sub-distribution family. The options
#'     are \code{"norm"} for FMR models, \code{"lnorm"} for mixture of AFT
#'     regression models with Log-Normal sub-distributions, \code{"weibull"}
#'     for mixture of AFT regression models with Weibull sub-distributions
#' @param initCoeff Vector of initial values for regression coefficients
#'     including intercepts
#' @param initDispersion Vector of initial values for standard deviations
#' @param initmixProp Vector of initial values for proportion of components
#' @param penFamily Penalty name that is used in variable selection method
#'     The available options are  \code{"lasso"}, \code{"adplasso"},
#'     \code{"mcp"}, \code{"scad"}, \code{"sica"} and \code{"hard"}.
#' @param lambPen A vector of positive numbers for tuning parameters
#' @param lambRidge A positive value for tuning parameter in Ridge
#'     Regression or Elastic Net
#' @param nIterEM Maximum number of iterations for EM algorithm
#' @param nIterNR Maximum number of iterations for Newton-Raphson algorithm
#' @param conveps A positive value for avoiding NaN in computing divisions
#' @param convepsEM A positive value for treshold of convergence in
#'     EM algorithm
#' @param convepsNR A positive value for treshold of convergence in
#'     NR algorithm
#' @param porNR A positive interger for maximum number of searches in
#' NR algorithm
#' @param gamMixPor Proportion of mixing parameters in the penalty. The
#'     value must be in the interval [0,1]. If \code{gamMixPor = 0}, the
#'     penalty structure is no longer mixture.
#' @param ... Other possible options
#' @keywords FMR AFT Censored EM Algorithm Ridge Regression
#' ElasticNet Selection LASSO MCP SCAD SICA Adaptive
#' @concept fmr, aft, lasso, adplasso, mcp, scad, sica, ridge, elastic net
#' @details The penalized likelihood of a finite mixture of AFT regression
#'     models is written as \deqn{\tilde\ell_{n}(\boldsymbol\Psi)
#'     =\ell_{n}(\boldsymbol\Psi) -
#'     \mathbf{p}_{\boldsymbol\lambda_{n}}(\boldsymbol\Psi)}
#'     where \deqn{\mathbf{p}_{\boldsymbol\lambda_{n}}(\boldsymbol\Psi) =
#'     \sum\limits_{k=1}^{K}\pi_{k}^\alpha\left\{
#'     \sum\limits_{j=1}^{d}p_{\lambda_{n,k}}(\beta_{kj}) \right\}.}
#'     In the M step of EM algorithm the
#'     function \deqn{\tilde{Q}(\boldsymbol\Psi,\boldsymbol\Psi^{(m)})
#'     =\sum\limits_{k=1}^{K} \tilde{Q}_{k}(\boldsymbol\Psi_k,
#'     \boldsymbol\Psi^{(m)}_k) =
#'     \sum\limits_{k=1}^{K} \left[{Q}_{k}(\boldsymbol\Psi_k,
#'     \boldsymbol\Psi^{(m)}_k) - \pi_{k}^\alpha\left\{
#'     \sum\limits_{j=1}^{d}p_{\lambda_{n,k}}(\beta_{kj}) \right\}\right]}
#'     is maximized. Since the penalty function is singular at origin, we
#'     use a local quadratic approximation (LQA) for the penalty as
#'     follows, \deqn{\mathbf{p}^\ast_{k,\boldsymbol\lambda_{n}}
#'     (\boldsymbol\beta,\boldsymbol\beta^{(m)})
#'     =(\pi_{k}^{(m)})^{\alpha}\sum\limits_{j=1}^{d}\left\{
#'     p_{\lambda_{n,k}}(\beta_{kj}^{(m)}) + { p^{\prime}_{\lambda_{n,k}}
#'     (\beta_{kj}^{(m)})  \over 2\beta_{kj}^{(m)}}(\beta_{kj}^{2} -
#'     {\beta_{kj}^{(m)}}^{2}) \right\}.} Therefore maximizing \eqn{Q} is
#'     equivalent to maximizing the
#'     function \deqn{ {Q}^\ast(\boldsymbol\Psi,\boldsymbol\Psi^{(m)})
#'     =\sum\limits_{k=1}^{K} {Q}^\ast_{k}(\boldsymbol\Psi_k,
#'     \boldsymbol\Psi^{(m)}_k) = \sum\limits_{k=1}^{K}
#'     \left[{Q}_{k}(\boldsymbol\Psi_k,\boldsymbol\Psi^{(m)}_k)-
#'     \mathbf{p}^\ast_{k,\boldsymbol\lambda_{n}}(\boldsymbol\beta,
#'     \boldsymbol\beta^{(m)})\right].}
#'     In case of Log-Normal sub-distributions, the maximizers of \eqn{Q_k}
#'     functions are as follows. Given the data and current estimates of
#'     parameters, the maximizers are \deqn{{\boldsymbol\beta}^{(m+1)}_{k}
#'     =({\boldsymbol z}^{\prime}\boldsymbol\tau^{(m)}_{k}{\boldsymbol z}+
#'     \varpi_{k}(\boldsymbol\beta_{kj}^{(m)}))^{-1}{\boldsymbol z}^{\prime}
#'     \boldsymbol\tau^{(m)}_{k}T^{(m)}_{k},}
#'     where \eqn{\varpi_{k}(\boldsymbol\beta_{kj}^{(m)})={diag}
#'     \left(\left(\pi_{k}^{(m+1)}\right)^\alpha
#'     \frac{{p}^{\prime}_{\lambda_{n},k}(\boldsymbol\beta_{kj}^{(m)})}
#'     {\boldsymbol\beta_{kj}^{(m)}}\right)}
#'     and \eqn{\sigma_{k}^{(m+1)}} is equal to \deqn{\sigma_{k}^{(m+1)}
#'     =\sqrt{\frac{\sum\limits_{i=1}^{n}\tau^{(m)}_{ik} (t^{(m)}_{ik}
#'     -{\boldsymbol z}_{i}\boldsymbol\beta^{(m)}_{k})^{2}}
#'     {\sum\limits_{i=1}^{n}\tau^{(m)}_{ik} {\left[\delta_{i}
#'     +(1-\delta_{i})\{A(w^{(m)}_{ik})[A(w^{(m)}_{ik})-
#'     w^{(m)}_{ik}]\}\right]}}}.}
#'     For the Weibull distribution, on the other hand,  we have
#'     \eqn{\tilde{\boldsymbol\Psi}^{(m+1)}_k
#'     =\tilde{\boldsymbol\Psi}^{(m)}_k
#'     - 0.5^{\kappa}\left[{H_{k}^{p,(m)}}\right]^{-1}I_{k}^{p,(m)}},
#'     where \eqn{H^p_{k}=H_k+h(\boldsymbol\Psi_k)}
#'     is the penalized version of hessian matrix
#'     and \eqn{I^p_{k}=I_k+h(\boldsymbol\Psi_k)\boldsymbol\Psi_k}
#'     is the penalized version of vector of first derivatives evaluated
#'     at \eqn{\tilde{\boldsymbol\Psi}_k^{(m)}}.
#' @references Shokoohi, F., Khalili, A., Asgharian, M. and Lin, S.
#' (2016 submitted) Variable Selection in Mixture of Survival Models
#' @return \code{\link{fmrsfit-class}}
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' nObs = 500
#' dispersion = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrs.gendata(nObs = nObs, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), dispersion = dispersion,
#'                      mixProp =mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#'
#' res.mle <- fmrs.mle(y = dat$y, x = dat$x, delta = dat$delta,
#'                    nComp = nComp, disFamily = "lnorm",
#'                    initCoeff = rnorm(nComp*nCov+nComp),
#'                    initDispersion = rep(1, nComp),
#'                    initmixProp = rep(1/nComp, nComp))
#'
#' res.lam <- fmrs.tunsel(y = dat$y, x = dat$x, delta = dat$delta,
#'                       nComp = ncomp(res.mle), disFamily = "lnorm",
#'                       initCoeff=c(coefficients(res.mle)),
#'                       initDispersion = dispersion(res.mle),
#'                       initmixProp = mixProp(res.mle),
#'                       penFamily = "adplasso")
#' res.var <- fmrs.varsel(y = dat$y, x = dat$x, delta = dat$delta,
#'                       nComp = ncomp(res.mle), disFamily = "lnorm",
#'                       initCoeff=c(coefficients(res.mle)),
#'                       initDispersion = dispersion(res.mle),
#'                       initmixProp = mixProp(res.mle),
#'                       penFamily = "adplasso",
#'                       lambPen = slot(res.lam, "lambPen"))
#'
#' coefficients(res.var)[-1,]
#' round(coefficients(res.var)[-1,],5)
#' @exportMethod fmrs.varsel
setGeneric("fmrs.varsel",
           function(y,delta,x,nComp,...) standardGeneric("fmrs.varsel"))


#' @title fmrs.gendata method
#' @description Generates a data set from Finite Mixture
#'     of AFT regression models or Finite Mixture of Regression models under
#'     the specified setting.
#' @author Farhad Shokoohi <shokoohi@icloud.com>
#' @family lnorm, norm, weibull
#' @name fmrs.gendata
#' @rdname fmrs.gendata-methods
#' @param nObs A numeric value represents number of observations (sample size)
#' @param nComp A numeric value represents the order (number of components)
#'     of an FMRs model
#' @param nCov A numberic value represents the number of covariates in
#'     design matrix
#' @param coeff A vector of all regression coefficients including
#'     intercepts. It must be a vector of length
#'     \code{nComp} by \code{nCov+1}.
#' @param dispersion A vector of positive values for dispersion parameters of
#'     sub-distributions in FMRs models
#' @param mixProp A vector of mixing proportions which their sum must be one
#' @param rho A numeric value in [-1, 1] which represents the correlation
#'     between covariates of design matrix
#' @param umax A numeric value represents the upper bound in Uniform
#'      distribution for censoring
#' @param disFamily A sub-distribution family. The options
#'     are \code{"lnormal"} for Log-Normal, \code{"norm"} for Normal and
#'     \code{"weibull"} for Weibull.
#' @param ... Other possible options
#' @import stats
#' @keywords FMRs AFT Censored Data Generation
#' @concept fmr, aft, censoring, data generation
#' @return A list including reponse, covariates and cenroing variables
#' @examples
#' set.seed(1980)
#' nComp = 2
#' nCov = 10
#' nObs = 500
#' REP = 500
#' dispersion = c(1, 1)
#' mixProp = c(0.4, 0.6)
#' rho = 0.5
#' coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
#' coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
#' umax = 40
#'
#' dat <- fmrs.gendata(nObs = nObs, nComp = nComp, nCov = nCov,
#'                      coeff = c(coeff1, coeff2), dispersion = dispersion,
#'                      mixProp =mixProp, rho = rho, umax = umax,
#'                      disFamily = "lnorm")
#' @exportMethod fmrs.gendata
setGeneric("fmrs.gendata",
           function(nObs,
                    nComp,
                    nCov,
                    coeff,
                    dispersion,
                    mixProp,
                    rho,
                    umax,...) standardGeneric("fmrs.gendata"))
