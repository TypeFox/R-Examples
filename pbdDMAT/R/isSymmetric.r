#' isSymmetric
#' 
#' Tests if a distributed matrix is symmetric.
#' 
#' The test is performed by comparing the object against its transpose.
#' 
#' @param object
#' Distributed matrix
#' @param tol
#' Numerical tolerance for the comparison.
#' @param ...
#' Additional arguments passed to \code{all.equal()}.
#' 
#' @name isSymmetric
#' @rdname isSymmetric
NULL

setGeneric(name = "isSymmetric", useAsDefault = base::isSymmetric, package="pbdDMAT")

#' @rdname isSymmetric
#' @export
setMethod("isSymmetric", signature(object="ddmatrix"), 
  function (object, tol = 100 * .Machine$double.eps, ...) 
  {
    if (object@dim[1L] != object@dim[2L]) 
      return(FALSE)
    
    test <- all.equal(object, t(object), tolerance = tol, ...)
    
    pbdMPI::comm.all(isTRUE(test))
  }
)





### TODO method for dsmatrix
