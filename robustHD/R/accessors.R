# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

## get component 'best'

getBest <- function(x, ...) UseMethod("getBest")

# get index of final model according to BIC
getBest.bicSelect <- function(x, ...) x$best

# get index of best reweighted and raw fit
getBest.fitSelect <- function(x, ...) x$best

# return NULL by default
getBest.default <- function(x, ...) NULL


## get a component for certain steps of the model sequence
# this is used for accessors that are exported to the namespace
# (coefficients, fitted values, residuals, ...), so checks for the
# arguments are necessary

getComponent <- function(x, which, ...) UseMethod("getComponent")

getComponent.seqModel <- function(x, component, s = NA,
                                  drop = !is.null(s), ...) {
  # extract component
  comp <- x[[component]]
  # check selected steps and extract corresponding parts of the component
  steps <- x$s  # computed steps
  if(length(steps) > 1) {
    if(!is.null(s)) {
      # check selected steps
      if(isTRUE(is.na(s))) s <- getSOpt(x)  # defaults to optimal step
      else s <- checkSteps(s, sMin=steps[1], sMax=steps[length(steps)], ...)
      # extract corresponding parts of the component
      if(is.null(dim(comp))) comp <- comp[s - steps[1] + 1]
      else comp <- comp[, s - steps[1] + 1, drop=FALSE]
    }
  }
  # drop dimension if requested and return component
  if(isTRUE(drop)) dropCol(comp) else comp
}

getComponent.sparseLTS <- function(x, which, s = NA,
                                   fit = c("reweighted", "raw", "both"),
                                   drop = !is.null(s), ...) {
  ## initializations
  fit <- match.arg(fit)
  if(fit != "reweighted") raw.which <- paste("raw", which, sep=".")
  sMax <- length(x$lambda)
  ## check lambda contains more than one value for the penalty parameter
  if(sMax > 1) {
    # extract component
    comp <- switch(fit, reweighted=x[[which]], raw=x[[raw.which]],
                   both={
                     rew <- x[[which]]
                     raw <- x[[raw.which]]
                     if(is.null(dim(rew))) unlist(list(reweighted=rew, raw=raw))
                     else {
                       colnames(rew) <- paste("reweighted",
                                              colnames(rew), sep=".")
                       colnames(raw) <- paste("raw", colnames(raw), sep=".")
                       cbind(rew, raw)
                     }
                   })
    # check selected steps and extract corresponding parts of the component
    if(!is.null(s)) {
      if(isTRUE(is.na(s))) {
        s <- getSOpt(x, fit=fit)  # defaults to optimal step
        if(fit == "both") s[2] <- sMax + s[2]
      } else if(fit == "both" && is.list(s)) {
        # list of steps for each fit
        s <- rep(s, length.out=2)
        s <- lapply(s, checkSteps, sMin=1, sMax=sMax)
        s <- c(s[[1]], sMax+s[[2]])
      } else {
        s <- checkSteps(s, sMin=1, sMax=sMax)
        if(fit == "both") s <- c(s, sMax+s)
      }
      # extract selected steps
      if(is.null(dim(comp))) comp <- comp[s]
      else comp <- comp[, s, drop=FALSE]
    }
    if(isTRUE(drop)) comp <- dropCol(comp)  # drop dimension if requested
  } else {
    # extract component
    comp <- switch(fit, reweighted=x[[which]], raw=x[[raw.which]],
                   both={
                     rew <- x[[which]]
                     raw <- x[[raw.which]]
                     cfun <- if(length(rew) > 1) cbind else c
                     cfun(reweighted=rew, raw=raw)
                   })
  }
  ## return component
  comp
}


#' Extract the residual scale of a robust regression model
#'
#' Extract the robust scale estimate of the residuals from a robust regression
#' model.
#'
#' Methods are implemented for models of class \code{"lmrob"} (see
#' \code{\link[robustbase]{lmrob}}), \code{"lts"} (see
#' \code{\link[robustbase]{ltsReg}}), \code{"rlm"} (see
#' \code{\link[MASS]{rlm}}), \code{"seqModel"} (see \code{\link{rlars}}) and
#' \code{"sparseLTS"} (see \code{\link{sparseLTS}}).  The default method
#' computes the MAD of the residuals.
#'
#' @param x  the model fit from which to extract the robust residual scale
#' estimate.
#' @param s  for the \code{"seqModel"} method, an integer vector giving
#' the steps of the submodels for which to extract the robust residual scale
#' estimate (the default is to use the optimal submodel).  For the
#' \code{"sparseLTS"} method, an integer vector giving the indices of the
#' models from which to extract the robust residual scale estimate.  If
#' \code{fit} is \code{"both"}, this can be a list with two components, with
#' the first component giving the indices of the reweighted fits and the second
#' the indices of the raw fits.  The default is to use the optimal model for
#' each of the requested estimators.  Note that the optimal models may not
#' correspond to the same value of the penalty parameter for the reweighted
#' and the raw estimator.
#' @param fit  a character string specifying from which fit to extract the
#' robust residual scale estimate.  Possible values are \code{"reweighted"}
#' (the default) for the residual scale of the reweighted fit, \code{"raw"} for
#' the residual scale of the raw fit, or \code{"both"} for the residual scale
#' of both fits.
#' @param \dots  additional arguments to be passed down to methods.
#'
#' @return
#' A numeric vector or matrix giving the robust residual scale estimates for
#' the requested model fits.
#'
#' @author Andreas Alfons
#'
#' @seealso \code{\link[=AIC.seqModel]{AIC}}, \code{\link[robustbase]{lmrob}},
#' \code{\link[robustbase]{ltsReg}}, \code{\link[MASS]{rlm}},
#' \code{\link{rlars}}, \code{\link{sparseLTS}}
#'
#' @examples
#' data("coleman")
#' fit <- lmrob(Y ~ ., data=coleman)
#' getScale(fit)
#'
#' @keywords regression
#'
#' @import stats
#' @export

getScale <- function(x, ...) UseMethod("getScale")

#' @method getScale default
#' @export
getScale.default <- function(x, ...) mad(residuals(x))  # use MAD by default

#' @method getScale lmrob
#' @export
getScale.lmrob <- function(x, ...) x$scale

#' @method getScale lts
#' @export
getScale.lts <- function(x, ...) x$scale

#' @method getScale rlm
#' @export
getScale.rlm <- function(x, ...) x$s

#' @rdname getScale
#' @method getScale seqModel
#' @export
getScale.seqModel <- function(x, s = NA, ...) getComponent(x, "scale", s=s, ...)

#' @rdname getScale
#' @method getScale sparseLTS
#' @export
getScale.sparseLTS <- function(x, s = NA,
                               fit = c("reweighted", "raw", "both"),
                               ...) {
  getComponent(x, "scale", s=s, fit=fit, ...)
}


## get optimal step

getSOpt <- function(x, ...) UseMethod("getSOpt")

getSOpt.seqModel <- function(x, ...) {
  sOpt <- getBest(x$crit)
  if(!is.null(sOpt)) sOpt <- sOpt + x$s[1] - 1
  sOpt
}

getSOpt.sparseLTS <- function(x, fit = "reweighted", ...) {
  sOpt <- getBest(x$crit)
  if(fit != "both") sOpt <- sOpt[fit]
  sOpt
}
