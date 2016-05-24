# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' Deprecated functions in package \pkg{robustHD}
#' 
#' These functions are provided for compatibility with older versions only, and 
#' may be defunct as soon as the next release.
#' 
#' \code{sparseLTSGrid} is a wrapper function for \code{\link{sparseLTS}} 
#' that only differs in the default values for the penalty parameter 
#' \code{lambda}.
#' 
#' @name robustHD-deprecated
#' 
#' @param formula  a formula describing the model.
#' @param data  an optional data frame, list or environment (or object coercible 
#' to a data frame by \code{\link{as.data.frame}}) containing the variables in 
#' the model.  If not found in data, the variables are taken from 
#' \code{environment(formula)}, typically the environment from which 
#' \code{sparseLTSGrid} is called.
#' @param x  a numeric matrix containing the predictor variables.
#' @param y  a numeric vector containing the response variable.
#' @param lambda  a numeric vector of non-negative values to be used as penalty 
#' parameter.
#' @param mode  a character string specifying the type of penalty parameter.  If 
#' \code{"lambda"}, \code{lambda} gives the grid of values for the penalty 
#' parameter directly.  If \code{"fraction"}, the smallest value of the penalty 
#' parameter that sets all coefficients to 0 is first estimated based on 
#' bivariate winsorization, then \code{lambda} gives the fractions of that 
#' estimate to be used (hence all values of \code{lambda} should be in the 
#' interval [0,1] in that case).
#' @param \dots  additional arguments to be passed down, eventually to 
#' \code{\link{sparseLTS}}.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link[base]{Deprecated}}
#' 
#' @keywords regression robust

NULL


#' @rdname robustHD-deprecated
#' @export

sparseLTSGrid <- function(x, ...) {
  .Deprecated("sparseLTS")
  UseMethod("sparseLTSGrid")
}


#' @rdname robustHD-deprecated
#' @method sparseLTSGrid formula
#' @export

sparseLTSGrid.formula <- function(formula, data, ...) {
  # get function call
  call <- match.call()
  call[[1]] <- as.name("sparseLTSGrid")
  # prepare model frame
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  if(is.empty.model(mt)) stop("empty model")
  # extract response and candidate predictors from model frame
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf)
  # check if the specified model contains an intercept
  # if so, remove the column for intercept and use 'intercept=TRUE'
  # otherwise use 'intercept=FALSE'
  whichIntercept <- match("(Intercept)", colnames(x), nomatch = 0)
  intercept <- whichIntercept > 0
  if(intercept) x <- x[, -whichIntercept, drop = FALSE]
  # call default method
  fit <- sparseLTSGrid(x, y, intercept=intercept, ...)
  fit$call <- call  # add call to return object
  fit$terms <- mt   # add model terms to return object
  fit
}


#' @rdname robustHD-deprecated
#' @method sparseLTSGrid default
#' @export

sparseLTSGrid.default <- function(x, y, lambda, mode = c("lambda", "fraction"), 
                                  ...) {
  # get function call
  call <- match.call()
  call[[1]] <- as.name("sparseLTSGrid")
  # initializations
  if(missing(lambda)) {
    # if penalty parameter is not supplied, use fractions of a robust estimate 
    # of the smallest value that sets all coefficients to zero
    lower <- if(nrow(x) > ncol(x)) 0 else 0.1
    lambda <- seq(from=1, to=lower, by=-0.1)
    mode <- "fraction"
  }
  fit <- sparseLTS(x, y, lambda=lambda, mode=mode, ...)
  fit$call <- call  # add call to return object
  fit
}
