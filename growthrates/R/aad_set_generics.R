#' Additional Generic Functions
#'
#' These functions are specifically defined for package \pkg{growthrates},
#' all other generics are imported.
#'
#' @param object name of a 'growthrate' object
#' @param \dots other arguments passed to the methods
#
#'
#' @rdname generics
#' @keywords internal
#' @include aac_classes.R
#'
#' @exportMethod rsquared
#'
setGeneric("rsquared", function(object, ...) standardGeneric("rsquared"))

#' @rdname generics
#' @exportMethod obs
#'
setGeneric("obs", function(object, ...) standardGeneric("obs"))

#' @rdname generics
#' @exportMethod results
#'
setGeneric("results", function(object, ...) standardGeneric("results"))

#' @rdname multisplit
#' @exportMethod multisplit
#'
setGeneric("multisplit", function(data, grouping, drop = TRUE, sep = ":", ...) standardGeneric("multisplit"))




