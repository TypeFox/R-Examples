#' LU Factorization
#' 
#' LU factorization for distributed matrices with R-like syntax, with
#' calculations performed by the PBLAS and ScaLAPACK libraries.
#' 
#' Extensions of R linear algebra functions.
#' 
#' @param x
#' numeric distributed matrices.
#' 
#' @return 
#' \code{lu()} performs LU factorization.
#' 
#' @examples
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:9, 3)
#' x <- as.ddmatrix(x)
#' 
#' y <- solve(t(A) %*% A)
#' print(y)
#' 
#' finalize()
#' }
#' 
#' @keywords Methods Linear Algebra
#' @aliases lu
#' @name ddmatrix-lu
#' @rdname ddmatrix-lu
NULL

setGeneric(name="lu", function(x) standardGeneric("lu"), package="pbdDMAT")

#' @rdname ddmatrix-lu
#' @export
setMethod("lu", signature(x="ddmatrix"), 
  function(x)
  {
    desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
    
    out <- base.rpdgetrf(a=x@Data, desca=desca)
    
    x@Data <- out
    
    return( x )
  }
)



