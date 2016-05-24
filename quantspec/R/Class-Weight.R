#' @include generics.R
#' @include Class-FreqRep.R
NULL

################################################################################
#' Interface Class to access different types of weighting functions.
#'
#' \code{Weights} is an S4 class that provides a common interface to
#' implementations of a weighting function \eqn{W_n(\omega)}{Wn(w)}.
#'
#' Currently three implementations are available:
#'     (1) \code{\link{KernelWeight}},
#'     (2) \code{\link{LagKernelWeight}} and
#'     (3) \code{\link{SpecDistrWeight}}.
#'
#' @name Weight-class
#' @aliases Weight
#'
#' @keywords S4-classes
#'
#' @slot values an array containing the weights.
#' @slot descr a description to be used in some plots.
#'
################################################################################
setClass(
    Class = "Weight",
    representation=representation(
        descr = "character"
    )
)

################################################################################
#' Get attribute \code{descr} from a \code{Weight}.
#'
#' @name getDescr-Weight
#' @aliases getDescr,Weight-method
#'
#' @keywords Access-functions
#'
#' @param object \code{Weight} from which to get the \code{descr}.
#' @return Returns the \code{descr} attribute.
################################################################################
setMethod(f = "getDescr",
    signature = "Weight",
    definition = function(object) {
      return(object@descr)
    }
)
