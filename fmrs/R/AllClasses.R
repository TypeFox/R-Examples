#' @title An S4 class to represent estimated optimal lambdas
#'
#' @description An S4 class to represent estimated optimal lambdas resulted
#'     from running \code{\link{fmrs.tunsel}}
#' @name fmrstunpar-class
#' @docType class
#' @exportClass fmrstunpar
#' @import methods
#' @slot ncomp A length-one numeric vector
#' @slot lambPen A dimension-one-\code{ncomp} numeric array
#' @slot disFamily A length-one character vector
#' @slot penFamily A length-one character vector
#' @slot lambRidge A length-one numeric vector
#' @slot model A length-one character vector
#' @rdname fmrstunpar-class
#' @keywords object
#' @export
fmrstunpar <- setClass("fmrstunpar",
                       representation(ncomp = "numeric",
                                      lambPen = "array",
                                      lambRidge = "numeric",
                                      disFamily = "character",
                                      penFamily = "character",
                                      model = "character"
                       ),
                       prototype (ncomp = numeric(length = 1L),
                                  lambPen = array(),
                                  lambRidge = numeric(length = 1L),
                                  disFamily = character(),
                                  penFamily = character(),
                                  model = character()
                       )
)

#' @title An S4 class to represent a fitted FMRs model
#'
#' @description is an S4 class represents a fitted of FMRs model
#'     resulted from running \code{\link{fmrs.mle}}
#'     or \code{\link{fmrs.varsel}}
#' @name fmrsfit-class
#' @import methods
#' @slot y A length-\code{nobs} numeric vector
#' @slot delta A length-\code{nobs} numeric vector
#' @slot x A dimension-\code{nobs}-\code{ncov} numeric matrix
#' @slot nobs A length-one numeric vector
#' @slot ncov A length-one numeric vector
#' @slot ncomp A length-one numeric vector
#' @slot coefficients A length-\code{(ncov+1)}-\code{ncomp} numeric matrix
#' @slot dispersion A length-\code{ncomp} numeric vector
#' @slot mixProp A length-\code{ncomp} numeric vector
#' @slot logLik A length-one numeric vector
#' @slot BIC A length-one numeric vector
#' @slot nIterEMconv A length-one numeric vector
#' @slot disFamily A length-one character vector
#' @slot penFamily A length-one character vector
#' @slot lambPen A length-\code{ncomp} numeric vector
#' @slot lambRidge A length-one numeric vector
#' @slot model A length-one character vector
#' @slot fitted A dimension-\code{nobs}-\code{ncomp} numeric matrix
#' @slot residuals A dimension-\code{nobs}-\code{ncomp} numeric matrix
#' @slot weights A dimension-\code{nobs}-\code{ncomp} numeric matrix
#' @docType class
#' @keywords object
#' @rdname fmrsfit-class
#' @exportClass fmrsfit
frmsfit <- setClass("fmrsfit",
                    representation(y = "vector",
                                   delta = "vector",
                                   x = "matrix",
                                   nobs = "numeric",
                                   ncov = "numeric",
                                   ncomp = "numeric",
                                   coefficients = "matrix",
                                   dispersion = "array",
                                   mixProp = "array",
                                   logLik = "numeric",
                                   BIC = "numeric",
                                   nIterEMconv = "numeric",
                                   disFamily = "character",
                                   penFamily = "character",
                                   lambPen = "array",
                                   lambRidge = "numeric",
                                   model = "character",
                                   fitted = "matrix",
                                   residuals = "matrix",
                                   weights = "matrix"
                    ),
                    prototype (y = vector(),
                               delta = vector(),
                               x = matrix(),
                               nobs = numeric(length = 1L),
                               ncov = numeric(length = 1L),
                               ncomp = numeric(length = 1L),
                               coefficients = matrix(),
                               dispersion = array(),
                               mixProp = array(),
                               logLik = numeric(length = 1L),
                               BIC = numeric(length = 1L),
                               nIterEMconv = numeric(length = 1L),
                               disFamily = character(),
                               penFamily = character(),
                               lambPen = array(),
                               lambRidge = numeric(length = 1L),
                               model = character(),
                               fitted = matrix(),
                               residuals = matrix(),
                               weights = matrix()
                    ),
                    contains = "fmrstunpar"
)



validity.fmrsfit <- function(object) {
  if( (length(object@y)!= object@nobs) | (dim(object@x)[1]!= object@nobs) |
      (dim(object@x)[1]!=length(object@y)) ) {
    return("data do not match")
  }else return(TRUE)
}


