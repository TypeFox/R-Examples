#' Row and Column binds for Distributed Matrices
#' 
#' Row and column binds.
#' 
#' The \code{...} list of arguments can be vectors, matrices, or distributed
#' matrices so long as non-distributed objects are not used with distributed
#' objects. This kind of mixing-and-matching will lead to chaos. Currently no
#' check is performed to prevent the user from this mixing-and-matching for
#' performance reasons (it is slow enough already).
#' 
#' @param ... 
#' vectors, matrices, or distributed matrices.
#' @param ICTXT 
#' BLACS communicator number for return object.
#' @param deparse.level 
#' integer controlling the construction of labels in the
#' case of non-matrix-like arguments. Does nothing for distributed matrices.
#' 
#' @return 
#' Returns a vector, matrix, or distributed matrix, depending on input.
#' 
#' @section Methods: \describe{ \item{list("signature(... = \"ANY\")")}{an R
#' object.} }
#' 
#' @examples
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' x <- ddmatrix(1:16, ncol=4, bldim=2)
#' 
#' y <- rbind(x, x)
#' 
#' print(y)
#' 
#' finalize()
#' }
#' 
#' @keywords Methods
#' @name binds
#' @rdname binds
NULL



### Bind 2 ddmatrix's
dmat.rbind2 <- function(args, ICTXT=.pbd_env$ICTXT)
{ 
  oldctxt <- args[[1]]@ICTXT
  
  args <- lapply(args, 
    FUN=function(dx) dmat.redistribute(dx=dx, bldim=dx@bldim, ICTXT=1)
  )
  
  dim <- c(sum(sapply(args, function(x) dim(x)[1])), args[[1]]@dim[2])
  bldim <- args[[1]]@bldim
  ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=1, fixme=TRUE)
  
  Data <- lapply(args, submatrix)
  
  ret <- new("ddmatrix", Data=Reduce(base::rbind, Data), dim=dim, ldim=ldim, bldim=bldim, ICTXT=1)
  
  if (ICTXT!=1)
    ret <- dmat.redistribute(dx=ret, bldim=ret@bldim, ICTXT=ICTXT)
  
  return( ret )
}



dmat.rbind <- function(..., ICTXT=.pbd_env$ICTXT)
{
  args <- list(...)
  
  return( dmat.rbind2(args=args, ICTXT=ICTXT) )
}



#' @rdname binds
#' @export
setMethod("rbind", "ANY", 
  function(..., ICTXT=.pbd_env$ICTXT, deparse.level=1)
  {
    args <- list(...)
    
    if (is.ddmatrix(args[[1]]))
      ret <- dmat.rbind2(args=args, ICTXT=ICTXT)
    else
      ret <- base::rbind(...=..., deparse.level=deparse.level)
    
    return( ret )
  }
)



dmat.cbind <- function(..., ICTXT=.pbd_env$ICTXT)
{
  args <- list(...)
  
  oldctxt <- args[[1]]@ICTXT
  
  args <- lapply(args, 
    FUN=function(dx) dmat.redistribute(dx=dx, bldim=dx@bldim, ICTXT=2)
  )
  
  dim <- c(args[[1]]@dim[1], sum(sapply(args, function(x) dim(x)[2])))
  bldim <- args[[1]]@bldim
  ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=2, fixme=TRUE)
  
  Data <- lapply(args, submatrix)
  
  ret <- new("ddmatrix", Data=Reduce(base::cbind, Data), dim=dim, ldim=ldim, bldim=bldim, ICTXT=2)
  
  if (ICTXT!=2)
    ret <- dmat.redistribute(dx=ret, bldim=ret@bldim, ICTXT=ICTXT)
  
  return( ret )
}



#' @rdname binds
#' @export
setMethod("cbind", "ANY", 
  function(..., ICTXT=.pbd_env$ICTXT, deparse.level=1)
  {
    args <- list(...)
    
    if (is.ddmatrix(args[[1]]))
      ret <- dmat.cbind(...=..., ICTXT=ICTXT)
    else
      ret <- base::cbind(...=..., deparse.level=deparse.level)
    
    return( ret )
  }
)
