# inversion via a cholesky, or inversion of crossprod(x) via qr

#' Inverse from Choleski (or QR) Decomposition
#' 
#' \code{qr()} takes the QR decomposition.
#' 
#' The function returns the inverse of a choleski factored matrix, or the
#' inverse of \code{crossprod(x)} if \code{qr.R(qr(x))} is passed.
#' 
#' @aliases chol2inv chol2inv-method chol2inv,ddmatrix-method chol2inv
#' @docType methods
#' @param x numeric distributed matrices for
#' @param size number of columns of \code{x} containing the Choleski
#' factorization.
#' @return A numeric distributed matrix.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{}
#' \item{list("signature(x = \"ANY\")")}{} }
#' @keywords Methods Linear Algebra
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' comm.set.seed(diff=T)
#' x <- ddmatrix("rnorm", 3, 3, bldim=2)
#' 
#' R <- qr.R(qr(x))
#' xtx.inv <- chol2inv(R)
#' 
#' id <- as.matrix(xtx.inv %*% crossprod(x))
#' 
#' comm.print(id)
#' 
#' finalize()
#' }
#' 
#' @name chol2inv
#' @rdname chol2inv
#' @export
setMethod("chol2inv", signature(x="ddmatrix"), 
  function(x, size = NCOL(x))
  {
    nr <- x@dim[1L]
    nc <- x@dim[2L]
    if (is.na(size) || size <= 0L || size > nr || size > nc) 
      comm.stop("invalid 'size' argument in 'chol2inv'")
    
    if (size < nr || size < nc){
      vec <- 1L:size
      x <- x[vec, vec]
    }
    
    descx <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
    
    cdim <- rep(size, 2)
    cldim <- base.numroc(dim=cdim, bldim=x@bldim, ICTXT=x@ICTXT)
    descc <- base.descinit(dim=cdim, bldim=x@bldim, cldim, ICTXT=x@ICTXT)
    
    out <- base.pdchtri(uplo='U', x=x@Data, descx=descx, descc=descc)
    
    c <- new("ddmatrix", Data=out, dim=cdim, ldim=cldim, bldim=x@bldim, ICTXT=x@ICTXT)
    
    return( c )
  }
)

