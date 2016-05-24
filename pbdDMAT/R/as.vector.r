#' Distributed object to Vector Converters
#' 
#' Converts a distributed matrix into a non-distributed vector.
#' 
#' The \code{proc.dest=} argument accepts either the BLACS grid position or the
#' MPI rank if the user desires a single process to own the matrix.
#' Alternatively, passing the default value of \code{'all'} will result in all
#' processes owning the matrix. If only a single process owns the undistributed
#' matrix, then all other processes store \code{NULL} for that object.
#' 
#' @param x 
#' numeric distributed matrix
#' @param ...
#' Additional arguments.
#' @param mode 
#' A character string giving an atomic mode or "list", or (except
#' for 'vector') "any".
#' @param proc.dest 
#' destination process for storing the matrix
#' 
#' @return Returns an ordinary R vector.
#' 
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' dx <- ddmatrix(1:16, ncol=4, bldim=2)
#' y <- as.vector(dx, proc.dest=0)
#' 
#' comm.print(y)
#' 
#' finalize()
#' }
#' 
#' @keywords Methods
#' @name as.vector
#' @rdname as.vector
setGeneric(name="as.vector", 
  function(x, ...)
    standardGeneric("as.vector"),
  useAsDefault=function(x, ...) 
    base::as.vector(x, ...)
)



#' @rdname as.vector
#' @export
setMethod("as.vector", signature(x="ddmatrix"), 
  function(x, mode="any", proc.dest="all"){
    ret <- as.vector(base.as.matrix(x, proc.dest=proc.dest), mode=mode)
    
    if (is.logical(x@Data))
      storage.mode(ret) <- "logical"
    
    return( ret )
  }
)
