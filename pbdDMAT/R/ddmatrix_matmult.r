#' Matrix Multiplication
#' 
#' Multiplies two distributed matrices, if they are conformable.
#' 
#' \code{x} and \code{y} must be conformable, on the same BLACS context, but
#' they need not be blocked with the same blocking dimension. The return will
#' default to the blocking dimension of \code{x}.
#' 
#' If you need to use \code{x} and \code{y} with differing blocking dimensions
#' and you want the return to have blocking different from that of \code{x},
#' then use the function \code{base.rpdgemm()}.
#' 
#' The \code{crossprod()} and \code{tcrossprod()} functions behave exactly as
#' their R counterparts.
#' 
#' @param x,y numeric distributed matrices
#' @return 
#' Returns a distributed matrix.
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
#' y <- x %*% x
#' print(y)
#' 
#' finalize()
#' }
#' 
#' @keywords Methods Linear Algebra
#' @name matmult
#' @rdname matmult
NULL

dmat.ddmatmult <- function(x, y, outbldim=x@bldim)
{
  if (!is.ddmatrix(x) || !is.ddmatrix(y))
    comm.stop("'x' and 'y' must be of class 'ddmatrix'")
  else if (x@dim[2L] != y@dim[1L])
    comm.stop("Error : non-conformable arguments.")
  
  base.checkem(x=x, y=y, checks=2)
  
  ICTXT <- x@ICTXT
  
  cdim <- c(x@dim[1L], y@dim[2L])
  cldim <- base.numroc(cdim, outbldim, ICTXT=ICTXT)
  
  descx <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=ICTXT)
  descy <- base.descinit(dim=y@dim, bldim=y@bldim, ldim=y@ldim, ICTXT=ICTXT)
  descc <- base.descinit(dim=cdim, bldim=outbldim, ldim=cldim, ICTXT=ICTXT)
  
  out <- base.rpdgemm(transx='N', transy='N', x=x@Data, descx=descx, y=y@Data, descy=descy, descc=descc)
  
  c <- new("ddmatrix", Data=out, dim=cdim, ldim=cldim, bldim=outbldim, ICTXT=ICTXT)
  
  return( c )
}

#' @rdname matmult
#' @export
setMethod("%*%", signature(x="ddmatrix", y="ddmatrix"),
  function(x, y)
    dmat.ddmatmult(x, y, outbldim=x@bldim)
)



dmat.crossprod <- function(trans, x)
{
  ICTXT <- x@ICTXT
  trans <- toupper(trans)
  
  if (trans=='N'){
    n <- x@dim[2L]
    k <- x@dim[1L]
  } else {
    n <- x@dim[1L]
    k <- x@dim[2L]
  }
  
  bldim <- x@bldim
  
  cdim <- c(n, n)
  cldim <- base.numroc(cdim, bldim, ICTXT=ICTXT)
  
  descx <- base.descinit(dim=x@dim, bldim=bldim, ldim=x@ldim, ICTXT=ICTXT)
  descc <- base.descinit(dim=cdim, bldim=bldim, ldim=cldim, ICTXT=ICTXT)
  
  out <- base.crossprod(uplo='U', trans=trans, x=x@Data, descx=descx, descc=descc)
  
  c <- new("ddmatrix", Data=out, dim=cdim, ldim=cldim, bldim=bldim, ICTXT=ICTXT)
  
  return( c )
}




#' @rdname matmult
#' @export
setMethod("crossprod", signature(x="ddmatrix", y="ANY"), 
  function(x, y = NULL)
  {
    if (is.null(y)){
      ret <- dmat.crossprod(trans='N', x=x)
      
      return( ret )
    }
    else if (!is.ddmatrix(y))
      comm.stop("Error : 'y' must be of class 'ddmatrix'.")
    else {
      if (x@dim[1L] != y@dim[1L])
        comm.stop("Error : non-conformable arguments.")
      
      base.checkem(x=x, y=y, checks=2)
      
      ret <- t(x) %*% y
      
      return( ret )
    }
  }
)



#' @rdname matmult
#' @export
setMethod("tcrossprod", signature(x="ddmatrix", y="ANY"), 
  function(x, y = NULL)
  {
    if (is.null(y)){
      ret <- dmat.crossprod(trans='T', x=x)
      
      return( ret )
    }
    else if (!is.ddmatrix(y))
      comm.stop("Error : 'y' must be of class 'ddmatrix'.")
    else {
      if (x@dim[1L] != y@dim[1L])
        comm.stop("Error : non-conformable arguments.")
      
      base.checkem(x=x, y=y, checks=2)
      
      ret <- x %*% t(y)
      
      return( ret )
    }
  }
)

