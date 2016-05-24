#####################################################
## Accessor methods                                ##
#####################################################

## Simple accessor for the (spherical) random effects.
##
## @title Random effects accessor
## @param object merMod object
## @param ... ignored
b <- function(object, ...) UseMethod("b")

## @rdname b
u <- function(object, ...) UseMethod("u")

#####################################################
## Summary / printing methods                      ##
#####################################################

##' The function \code{getInfo} is internally used to 
##' prepare object for producing a comparison chart in
##' \code{compare}.
##'
##' @rdname compare
##' @param object object
##' @return \code{getInfo} returns alist with estimated coefficients,
##' estimated variance components, sigma, deviance and parameter
##' configuration used to fit.
##' @examples
##' str(getInfo(fm1))
##' @export
getInfo <- function(object, ...) UseMethod("getInfo")
