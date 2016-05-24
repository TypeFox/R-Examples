## -------------------------------------------------------------
## methods for top-level class growthrates_fit
## -------------------------------------------------------------


#' Accessor Methods of Package \pkg{growthrates}.
#'
#' Functions to access the results of fitted growthrate objects:  \code{summary},
#'  \code{coeff}, \code{rsquared}, \code{deviance}, \code{residuals},
#'  \code{df.residual}, \code{obs}, \code{results}.
#'
#' @param object name of a 'growthrate' object.
#' @param cov boolean if the covariance matrix should be printed.
#' @param \dots other arguments passed to the methods.
#'
#' @rdname methods
#' @export rsquared
#' @exportMethod rsquared
#'
#' @examples
#'
#' data(bactgrowth)
#' splitted.data <- multisplit(bactgrowth, c("strain", "conc", "replicate"))
#'
#' ## get table from single experiment
#' dat <- splitted.data[[10]]
#'
#' fit1 <- fit_spline(dat$time, dat$value, spar=0.5)
#' coef(fit1)
#' summary(fit1)
#'
#' ## derive start parameters from spline fit
#' p <- c(coef(fit1), K = max(dat$value))
#' fit2 <- fit_growthmodel(grow_logistic, p=p, time=dat$time, y=dat$value, transform="log")
#' coef(fit2)
#' rsquared(fit2)
#' deviance(fit2)
#'
#' summary(fit2)
#'
#' plot(residuals(fit2) ~ obs(fit2)[,2])
#'
#'
setMethod("rsquared", "growthrates_fit",
          function(object, ...) {
            object@rsquared
          }
)


#' @rdname methods
#' @exportMethod obs
#'
setMethod("obs", "growthrates_fit",
          function(object, ...) {
            object@obs
          }
)


## -------------------------------------------------------------
## methods for sub-class nonlinear_fit
## -------------------------------------------------------------

#' @rdname methods
#' @exportMethod coef
#'
setMethod("coef", "growthrates_fit",
          function(object, ...) {
            #coef(object@fit, ...)
            object@par
          }
)

## necessary because intercept in fit (=lm) is log-transformed
#' @rdname methods
#' @exportMethod coef
#'
setMethod("coef", "easylinear_fit",
          function(object, ...) {
            object@par
          }
)


#' @rdname methods
#' @exportMethod coef
#'
setMethod("coef", "smooth.spline_fit",
          function(object, ...) {
            object@par
          }
)

#' @rdname methods
#' @exportMethod deviance
#'
setMethod("deviance", "growthrates_fit",
          function(object, ...) {
            deviance(object@fit, ...)
          }
)

#' @rdname methods
#' @exportMethod summary
#'
setMethod("summary", "growthrates_fit",
          function(object, ...) {
            summary(object@fit, cov=cov, ...)
          }
)

#' @rdname methods
#' @exportMethod summary
#'
setMethod("summary", "nonlinear_fit",
          function(object, cov=TRUE, ...) {
            summary(object@fit, cov=cov, ...)
          }
)


#' @rdname methods
#' @exportMethod residuals
#'
setMethod("residuals", "growthrates_fit",
          function(object, ...) {
            residuals(object@fit, ...)
          }
)

#' @rdname methods
#' @exportMethod df.residual
#'
setMethod("df.residual", "growthrates_fit",
          function(object, ...) {
            df.residual(object@fit, ...)
          }
)

## -------------------------------------------------------------
## methods for sub-class smooth.spline_fit
## -------------------------------------------------------------

#' @rdname methods
#' @exportMethod summary
#'
setMethod("summary", "smooth.spline_fit",
          function(object, cov=TRUE, ...) {
            coef <- coef(object)
            xy   <- object@xy
            cat("Fitted smoothing spline:\n")
            print(object@fit, ...)
            cat("\n")
            cat("Parameter values of exponential growth curve:\n")

            cat("Maximum growth at x=", xy[1], ", y=", xy[2], "\n")
            cat("y0 =", coef["y0"], "\n")
            cat("mumax =", coef["mumax"], "\n")
            cat("\n")
            cat("r2 of log transformed data=", rsquared(object), "\n")
          }
)


#' @rdname methods
#' @exportMethod df.residual
#'
setMethod("df.residual", "smooth.spline_fit",
          function(object, ...) {
            object@fit$df
          }
)

#' @rdname methods
#' @exportMethod deviance
#'
setMethod("deviance", "smooth.spline_fit",
          function(object, ...) {
            sum(residuals(object)^2)
          }
)


### ============================================================================
### Methods for Multiple Fits
### ============================================================================

#' @rdname methods
#' @exportMethod coef
#'
setMethod("coef", "multiple_fits",
          function(object, ...) {
            t(sapply(object@fits, function(x) coef(x, ...)))
          })

#' @rdname methods
#' @exportMethod rsquared
#'
setMethod("rsquared", "multiple_fits",
          function(object, ...) {
            sapply(object@fits, function(x) rsquared(x, ...))
          })


#' @rdname methods
#' @exportMethod deviance
#'
setMethod("deviance", "multiple_fits",
          function(object, ...) {
            sapply(object@fits, function(x) deviance(x, ...))
          })

#' @rdname methods
#' @exportMethod results
#'
setMethod("results", "multiple_fits",
          function(object, ...) {
            grouping <- object@grouping
            ret <- cbind(coef(object, ...), r2=rsquared(object, ...))
            keys <- matrix(unlist(strsplit(row.names(ret), ":")),
                           ncol=length(grouping), byrow=TRUE)

            keys <- as.data.frame(keys, stringsAsFactors = FALSE)

            ## try to convert keys to numeric or factor
            keys <- as.data.frame(lapply(keys, type.convert))

            ret <- cbind(keys, ret)
            names(ret)[1:length(grouping)] <- c(grouping)
            ret
          })


#' @rdname methods
#' @exportMethod results
#'
setMethod("results", "multiple_easylinear_fits",
          function(object, ...) {
            ret <- cbind(coef(object, ...), r2=rsquared(object, ...))
            keys <- matrix(unlist(strsplit(row.names(ret), ":")),
                           ncol=length(object@grouping), byrow=TRUE)

            keys <- as.data.frame(keys, stringsAsFactors = FALSE)

            ## try to convert keys to numeric or factor
            keys <- as.data.frame(lapply(keys, type.convert))

            ret <- cbind(keys, ret)
            names(ret) <- c(object@grouping, "y0", "mumax", "r2")
            ret
          })


#' @rdname methods
#' @exportMethod summary
#'
setMethod("summary", "multiple_fits",
          function(object, ...) {
            lapply(object@fits, summary)
          })


#' @rdname methods
#' @exportMethod residuals
#'
setMethod("residuals", "multiple_fits",
          function(object, ...) {
            unlist(lapply(object@fits, residuals))
          })

