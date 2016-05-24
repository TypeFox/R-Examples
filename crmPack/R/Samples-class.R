#####################################################################################
## Author: Daniel Sabanes Bove [sabanesd *a*t* roche *.* com]
## Project: Object-oriented implementation of CRM designs
##
## Time-stamp: <[Samples-class.R] by DSB Fre 16/01/2015 12:05>
##
## Description:
## A class for the MCMC samples. We need to have something slightly more
## flexible than the "mcmc" class from the "coda" package
##
## History:
## 25/03/2014   file creation
#####################################################################################

##' @include McmcOptions-class.R
##' @include McmcOptions-methods.R
{}

## --------------------------------------------------
## Class for the MCMC output
## --------------------------------------------------

##' Class for the MCMC output
##'
##' @slot data a list where each entry contains the samples of a (vector-valued)
##' parameter in a vector/matrix in the format (number of samples) x (dimension
##' of the parameter).
##' @slot options the \code{\linkS4class{McmcOptions}} which have been used
##'
##' @example examples/Sample-class-Samples.R
##' @export
##' @keywords classes
.Samples <-
    setClass(Class="Samples",
             representation(data="list",
                            options="McmcOptions"),
             prototype(data=
                           list(alpha=matrix(0, nrow=1, ncol=1),
                                beta=matrix(0, nrow=1, ncol=1)),
                       options=
                           McmcOptions(burnin=1,
                                       step=1,
                                       samples=1)),
             validity=
                 function(object){
                     o <- Validate()

                     o$check(all(sapply(object@data,
                                        NROW) == sampleSize(object@options)),
                             "all data elements must have as many rows as the sample size was")

                     o$result()
                 })
validObject(.Samples())


##' Initialization function for "Samples"
##'
##' @param data see \code{\linkS4class{Samples}}
##' @param options see \code{\linkS4class{Samples}}
##' @return the \code{\linkS4class{Samples}} object
##'
##' @export
##' @keywords methods
Samples <- function(data,
                    options)
{
    .Samples(data=data,
             options=options)
}

