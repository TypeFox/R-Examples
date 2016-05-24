#' Solve
#' 
#' Solving linear systems and matrix inversion for distributed matrices with R-like syntax, with
#' calculations performed by the PBLAS and ScaLAPACK libraries.
#' 
#' Extensions of R linear algebra functions.
#' 
#' @param a,b
#' numeric distributed matrices.  Here, \code{a}
#' and \code{b} must be on the same BLACS context and have the same blocking
#' dimension.
#' 
#' @return 
#' \code{solve()} solves systems and performs matrix inversion when argument
#' \code{b=} is missing.
#' 
#' @examples
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' x <- ddmatrix(1:9, 3)
#' 
#' y <- solve(t(A) %*% A)
#' print(y)
#' 
#' finalize()
#' }
#' 
#' @keywords Methods Linear Algebra
#' @name ddmatrix-solve
#' @rdname ddmatrix-solve
NULL



### inversion
#' @rdname ddmatrix-solve
#' @export
setMethod("solve", signature(a="ddmatrix"), 
  function(a)
  {
    if (diff(a@dim)!=0)
      comm.stop(paste("'a' (", a@dim[1], " x ", a@dim[2], ") must be square", sep=""))
    
    desca <- base.descinit(dim=a@dim, bldim=a@bldim, ldim=a@ldim, ICTXT=a@ICTXT)
    
    n <- desca[4L]
    
    out <- base.rpdgetri(n=n, a=a@Data, desca=desca)
    
    a@Data <- out
    
    return(a)
  }
)

### Solving systems
#' @rdname ddmatrix-solve
#' @export
setMethod("solve", signature(a="ddmatrix", b="ddmatrix"), 
  function(a, b)
  {
    base.checkem(x=a, y=b, checks=2:3)
    
    # Matrix descriptors
    desca <- base.descinit(dim=a@dim, bldim=a@bldim, ldim=a@ldim, ICTXT=a@ICTXT)
    descb <- base.descinit(dim=b@dim, bldim=b@bldim, ldim=b@ldim, ICTXT=b@ICTXT)
    
    n <- desca[4L]
    nrhs <- descb[4L]
    
    ret <- base.rpdgesv(n=n, nrhs=nrhs, a=a@Data, desca=desca, b=b@Data, descb=descb)
    
    b@Data <- ret
    
    return( b )
  }
)


