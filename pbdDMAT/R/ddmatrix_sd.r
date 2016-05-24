#' Covariance and Correlation
#' 
#' \code{sd()} will compute the standard deviations of the columns, equivalent
#' to calling \code{apply(x, MARGIN=2, FUN=sd)} (which will work for
#' distributed matrices, by the way). However, this should be much faster and
#' use less memory than \code{apply()}.  If \code{reduce=FALSE} then the return
#' is a distributed matrix consisting of one (global) row; otherwise, an
#' \code{R} vector is returned, with ownership of this vector determined by
#' \code{proc.dest}.
#' 
#' @param x
#' numeric distributed matrices.
#' @param na.rm
#' Logical; if TRUE, then \code{na.exclude()} is called first.
#' @param reduce 
#' logical or string. See details
#' @param proc.dest 
#' Destination process (or 'all') if a reduction occurs
#' 
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
#' x <- ddmatrix("rnorm", nrow=3, ncol=3)
#' 
#' cv <- cov(x)
#' print(cv)
#' 
#' finalize()
#' }
#' 
#' @keywords Methods
#' @name sd
#' @rdname sd
NULL



setGeneric(name="sd", 
  function(x, ...)
    standardGeneric("sd"),
  package="pbdDMAT"
)


#' @rdname sd
#' @export
setMethod("sd", signature(x="ddmatrix"),
function (x, na.rm = FALSE, reduce = FALSE, proc.dest="all") 
  {
    if (na.rm)
      x <- na.exclude(x)
    
    descx <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
    sdv <- base.pdclvar(x=x@Data, descx=descx)
    
    sdv <- sqrt(sdv)
    dim(sdv) <- c(1L, base::length(sdv))
#    sdv <- matrix(sqrt(sdv), nrow=1)
    
    ret <- new("ddmatrix", Data=sdv, dim=c(1, x@dim[2]), ldim=dim(sdv), bldim=x@bldim, ICTXT=x@ICTXT)
    
    if (reduce)
      ret <- as.vector(x=ret, proc.dest=proc.dest)
    
    return( ret )
  }
)

#' @rdname sd
#' @export
setMethod("sd", signature(x="ANY"), 
  function(x, na.rm = FALSE) 
    stats::sd(x=x, na.rm=na.rm)
)



