#' The ecd class
#' 
#' This S4 class is the major object class for elliptic distribution.
#' It stores the ecd parameters, numerical constants that facilitates
#' quadpack integration, statistical attributes, and optionally,
#' an internal structure for the quantile function.
#'
#' @name ecd-class
#'
#' @slot call The match.call slot
#' @slot alpha,gamma,sigma,beta,mu a length-one numeric.
#'                                 These are core ecd parameters.
#' @slot cusp a length-one numeric as cusp indicator. 
#'            0: not a cusp; 
#'            1: cusp specified by \code{alpha}; 
#'            2: cusp specified by \code{gamma}.
#' @slot lambda a length-one numeric, the leading exponent for the special model, default is 3.
#' @slot R,theta a length-one numeric. These are derived ecd parameters in polar coordinate.
#' @slot use.mpfr logical, internal flag indicating whether to use mpfr.
#' @slot const A length-one numeric as the integral of \eqn{exp(y(x))} that normalizes the PDF.
#' @slot const_left_x A length-one numeric marking the left point of PDF integration.
#' @slot const_right_x A length-one numeric marking the right point of PDF integration.
#' @slot stats A list of statistics, see \code{ecd.stats} for more information.
#' @slot quantile  An object of ecdq class, for quantile calculation.
#' @slot model A vector of four strings representing internal classification:
#'             \code{long_name.skew}, code{long_name},
#'             \code{short_name.skew}, \code{short_name}.
#'             This slot doesn't have formal use yet.
#'
#' @keywords class constructor
#'
#' @include ecd-package.R
#' @include ecd-numericMpfr-class.R
#' @include ecdq-class.R
#'
#' @exportClass ecd
setClass("ecd",
         representation(call = "call",
                        alpha = "numericMpfr",
                        gamma = "numericMpfr",
                        sigma = "numericMpfr",
                        beta  = "numericMpfr",
                        mu    = "numericMpfr",
                        cusp = "numericMpfr",
                        lambda = "numericMpfr",
                        R = "numericMpfr",
                        theta = "numericMpfr",
                        use.mpfr = "logical",
                        const = "numericMpfr",
                        const_left_x = "numericMpfr",
                        const_right_x = "numericMpfr",
                        stats = "list",
                        quantile = "ecdq",
                        model = "character"),
          prototype(call = call("ecd"),
                    alpha = 0,
                    gamma = 0,
                    sigma = 1,
                    beta  = 0,
                    mu    = 0,
                    cusp  = 0,
                    lambda = 3,
                    R = 0,
                    theta = 0,
                    use.mpfr = FALSE,
                    const = NaN,
                    const_left_x = NaN,
                    const_right_x = NaN,
                    stats = list(),
                    quantile = NULL,
                    model = "Elliptic")
)
### <---------------------------------------------------------------------->
setGeneric("const", function(object) standardGeneric("const"))

setMethod("const", "ecd", function(object){ object@const })

setGeneric("const<-", function(object, value) standardGeneric("const<-"))

setMethod("const<-", "ecd", function(object, value){
    
    object@const <- unname(value$const)
    object@const_left_x <- unname(value$const_left_x)
    object@const_right_x <- unname(value$const_right_x)
    
    stopifnot(is.numericMpfr(object@const))
    stopifnot(is.numericMpfr(object@const_left_x))
    stopifnot(is.numericMpfr(object@const_right_x))
    
    object
})
### <---------------------------------------------------------------------->
# don't set getter for stats, it interferes with stats package

setGeneric("stats<-", function(object, value) standardGeneric("stats<-"))

setMethod("stats<-", "ecd", function(object, value){
    object@stats <- value
    invisible(object)
})
### <---------------------------------------------------------------------->
# The setter for quantile data object ecdq

setGeneric("quantile<-", function(object, value) standardGeneric("quantile<-"))

setMethod("quantile<-", "ecd", function(object, value){
    object@quantile <- value
    invisible(object)
})
### <---------------------------------------------------------------------->
# cusp: 0 - not a cusp; 
#       1 - cusp with alpha specified; 2 - cusp with gamma specified

setGeneric("cusp<-", function(object, value) standardGeneric("cusp<-"))

setMethod("cusp<-", "ecd", function(object, value){
    object@cusp <- value
    object
})
### <---------------------------------------------------------------------->
