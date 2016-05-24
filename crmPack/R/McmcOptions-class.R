#####################################################################################
## Author: Daniel Sabanes Bove [sabanesd *a*t* roche *.* com]
## Project: Object-oriented implementation of CRM designs
##
## Time-stamp: <[McmcOptions-class.R] by DSB Fre 16/01/2015 12:02>
##
## Description:
## Encapsulate the three canonical MCMC options in a formal class.
##
## History:
## 12/12/2013   file creation (copy from glmBfp package)
## 19/12/2013   simplify doc (roxygen2 3.0)
## 29/01/2014   use initialize method rather than constructor function
##              to be better extensible
###################################################################################

##' @include helpers.R
{}


##' Class for the three canonical MCMC options
##'
##' @slot iterations number of MCMC iterations
##' @slot burnin number of burn-in iterations which are not saved
##' @slot step only every step-th iteration is saved after the burn-in
##'
##' @example examples/McmcOptions-class-McmcOptions.R
##' @export
##' @keywords classes
.McmcOptions <-
    setClass(Class="McmcOptions",
             representation(iterations="integer",
                            burnin="integer",
                            step="integer"),
             prototype(iterations=1000L,
                       burnin=100L,
                       step=2L),
             validity=function(object){
                 o <- Validate()

                 o$check(is.scalar(object@burnin) && (object@burnin >= 0L),
                         "burn-in must be non-negative scalar")
                 o$check(is.scalar(object@iterations),
                         "iterations must be integer scalar")
                 o$check(object@burnin < object@iterations,
                         "burn-in must be smaller than iterations")
                 o$check(is.scalar(object@step) && (object@step >= 1),
                         "step size must be scalar of at least 1")

                 o$result()
             })
validObject(.McmcOptions())


##' Initialization function for the "McmcOptions" class
##'
##' @param burnin number of burn-in iterations which are not saved (default:
##' \code{10,000})
##' @param step only every step-th iteration is saved after the burn-in
##' (default: \code{2})
##' @param samples number of resulting samples (by default \code{10,000} will
##' result)
##' @return the \code{\linkS4class{McmcOptions}} object
##'
##' @export
##' @keywords methods
McmcOptions <- function(burnin=1e4L,
                        step=2L,
                        samples=1e4L)
{
    .McmcOptions(iterations=safeInteger(burnin + (step * samples)),
                 burnin=safeInteger(burnin),
                 step=safeInteger(step))
}
validObject(McmcOptions())
