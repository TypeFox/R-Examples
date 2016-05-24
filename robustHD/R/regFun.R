# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

## get a call function
# this returns a function that either
# 1) simply evaluates a supplied function for the basic arguments if there are
#    no additional arguments in list format
# 2) evaluates a supplied function with 'do.call' if there are additional 
#    arguments in list format
getCallFun <- function(args) {
  if(length(args) == 0) function(..., fun, args) fun(...)
  else function(..., fun, args) do.call(fun, c(list(...), args))
}

## get the control object for model functions
#' @import robustbase MASS
getRegControl <- function(fun) {
  if(identical(fun, lmrob)) {
    fun <- .lmrob.fit
    useFormula <- FALSE
  } else if(identical(fun, ltsReg)) {
    fun <- .ltsReg
    useFormula <- FALSE
  } else if(identical(fun, rlm)) {
    fun <- .rlm
    useFormula <- FALSE
  } else if(identical(fun, lm)) {
    fun <- .lm.fit
    useFormula <- FALSE
  } else useFormula <- TRUE
  list(fun=fun, useFormula=useFormula)
}

## wrapper function for lm.fit() that returns appropriate class
.lm.fit <- function(x, y, ...) {
  fit <- lm.fit(x, y, ...)
  class(fit) <- "lm"
  fit
}

## wrapper function for lmrob.fit() with different default arguments
.lmrob.fit <- function(x, y, control, max.it = 500, k.max = 2500, ...) {
  if(missing(control)) {
    control <- lmrob.control(max.it=max.it, k.max=k.max, ...)
  }
  lmrob.fit(x, y, control=control)
}

## wrapper function for ltsReg() that handles constant column for intercept
.ltsReg <- function(x, y, intercept, ...) {
  ltsReg(removeIntercept(x), y, intercept=TRUE, ...)
}

## wrapper function for rlm() with different default arguments
.rlm <- function(x, y, maxit = 500, ...) rlm(x, y, maxit=maxit, ...)
