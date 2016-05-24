#####################################################################################
## Author: Daniel Sabanes Bove [sabanesd *a*t* roche *.* com]
## Project: Object-oriented implementation of CRM designs
##
## Time-stamp: <[McmcOptions-methods.R] by DSB Fre 16/01/2015 11:58>
##
## Description:
## Functions/methods for the MCMC formal class.
##
## History:
## 12/12/2013   file creation (copy from glmBfp package)
#####################################################################################

##' Compute the number of samples for a given MCMC options triple
##'
##' @param mcmcOptions the \code{\linkS4class{McmcOptions}} object
##' @return the resulting sample size
##'
##' @example examples/McmcOptions-methods-sampleSize.R
##' @export
##' @keywords programming
sampleSize <-
    function(mcmcOptions)
{
    stopifnot(is(mcmcOptions, "McmcOptions"))

    return(as.integer(ceiling((mcmcOptions@iterations - mcmcOptions@burnin) /
                                  mcmcOptions@step)))
}

##' Determine if we should save this sample
##'
##' @param iteration the current iteration index
##' @param mcmcOptions  the \code{\linkS4class{McmcOptions}} object
##' @return Logical value, if we should save this sample
##'
##' @example examples/McmcOptions-methods-saveSample.R
##' @export
##' @keywords programming internal
saveSample <-
    function(iteration,
             mcmcOptions)
{
    stopifnot(is(mcmcOptions, "McmcOptions"))

    return((iteration > mcmcOptions@burnin) &&
           (((iteration - mcmcOptions@burnin) %% mcmcOptions@step) == 0))
}



