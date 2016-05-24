##' The comparison operators \code{>}, \code{<}, \code{>=}, \code{<=},
##' \code{==}, and \code{!=} for \code{hyperSpec} objects.
##' 
##' \code{all.equal} checks the equality of two hyperSpec objects.
##' 
##' The comparison operators \code{>}, \code{<}, \code{>=}, \code{<=},
##' \code{==}, and \code{!=} work on the spectra matrix of the \code{hyperSpec}
##' object. They have their usual meaning (see \code{\link[base]{Comparison}}).
##' The operators work also with one \code{hyperSpec} object and a numeric
##' (scalar) object or a matrices of the same size as the spectra matrix of the
##' \code{hyperSpec} object.
##' 
##' With numeric vectors \code{\link[hyperSpec]{sweep}} might be more
##' appropriate.
##' 
##' If you want to calculate on the \code{data.frame} \code{hyperSpec@@data},
##' you have to do this directly on \code{hyperSpec@@data}.
##'
##' @author C. Beleites
##' @title Comparison of hyperSpec objects
##' @name Comparison
##' @rdname Comparison
##' @docType methods
##' @aliases Comparison Operators Compare,hyperSpec-method
##' Compare,hyperSpec,hyperSpec-method <,hyperSpec,hyperSpec-method >,hyperSpec,hyperSpec-method
##' <=,hyperSpec,hyperSpec-method >=,hyperSpec,hyperSpec-method ==,hyperSpec,hyperSpec-method
##' !=,hyperSpec,hyperSpec-method Compare,hyperSpec,matrix-method Compare,hyperSpec,numeric-method
##' Compare,matrix,hyperSpec-method Compare,numeric,hyperSpec-method
##' @param e1,e2 Either two \code{hyperSpec} objects or one \code{hyperSpec}
##'   object and matrix of same size as \code{hyperSpec[[]]} or a scalar
##'   (numeric of length 1).
##' 
##' As \code{hyperSpec} objects must have numeric spectra matrices, the
##'   resulting matrix of the comparison is returned directly.
##' @return a logical matrix for the comparison operators.
##' @seealso \code{\link[hyperSpec]{sweep-methods}} for calculations involving
##'   a vector and the spectral matrix.
##' 
##' \code{\link[methods]{S4groupGeneric}} for group generic methods.
##' 
##' \code{\link[base]{Comparison}} for the base comparison functions.
##' 
##' \code{\link[hyperSpec]{Arith}} for arithmetic operators,
##'   \code{\link[hyperSpec]{Math}} for mathematical group generic functions
##'   (groups Math and Math2) working on \code{hyperSpec} objects.
##' @keywords methods arith
##' @export
##' @examples
##' 
##' flu [,,445 ~ 450] > 300
##' 
##' all (flu == flu[[]])
##' 

setMethod ("Compare", signature (e1 = "hyperSpec", e2 = "hyperSpec"),
           function (e1, e2){
             validObject (e1)
             validObject (e2)

             callGeneric (e1[[]], e2[[]])
           }
           )

.compx <- function (e1, e2){
  validObject (e1)
  callGeneric (e1 [[]], e2)
}

.compy <- function (e1, e2){
  validObject (e2)
  callGeneric (e1, e2 [[]])
}

##' @rdname Comparison
setMethod ("Compare", signature (e1 = "hyperSpec", e2 = "numeric"), .compx)
##' @rdname Comparison
setMethod ("Compare", signature (e1 = "hyperSpec", e2 = "matrix"), .compx)

##' @rdname Comparison
setMethod ("Compare", signature (e1 = "numeric", e2 = "hyperSpec"), .compy)
##' @rdname Comparison
setMethod ("Compare", signature (e1 = "matrix", e2 = "hyperSpec"), .compy)


