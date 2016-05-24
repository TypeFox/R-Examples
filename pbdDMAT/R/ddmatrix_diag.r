#' Distributed Matrix Diagonals
#' 
#' Get the diagonal of a distributed matrix, or construct a distributed matrix
#' which is diagonal.
#' 
#' Gets the diagonal of a distributed matrix and stores it as a global R vector
#' owned by all processes.
#' 
#' @param x 
#' distributed matrix or a vector.
#' @param nrow,ncol 
#' in the case that \code{x} is a vector, these specify the
#' global dimension of the diagonal distributed matrix to be created.
#' @param type 
#' character. Options are 'matrix' or 'ddmatrix', with partial
#' matching.  This specifies the return type.
#' @param ... 
#' Extra arguments
#' @param min,max 
#' Min and max values for random uniform generation.
#' @param mean,sd 
#' Mean and standard deviation for random normal generation.
#' @param rate 
#' Rate for random exponential generation.
#' @param shape,scale 
#' Shape and scale parameters for random weibull generation.
#' @param bldim 
#' blocking dimension.
#' @param ICTXT 
#' BLACS context number.
#' 
#' @return 
#' If a distributed matrix is passed to \code{diag()} then it returns a
#' global R vector.
#' 
#' If a vector (numeric or character) is passed to \code{diag()} and
#' \code{type='ddmatrix'}, then the return is a diagonal distributed matrix.
#' 
#' @examples
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' ### Grab diagonal of a ddmatrix
#' x <- ddmatrix(1:16, 4, bldim=2)
#' y <- diag(x)
#' comm.print(y)
#' 
#' ### Construct a diagonal ddmatrix
#' eye <- diag(1, 10, type="ddmatrix", bldim=2)
#' eye
#' 
#' finalize()
#' }
#' 
#' @keywords Methods Extraction
#' @aliases diag
#' @name diag-constructors
#' @rdname diag-constructors
NULL



setGeneric(name="diag",
  function(x, ...)
    standardGeneric("diag"),
  package="pbdDMAT"
)



#' @rdname diag-constructors
#' @export
setMethod("diag", signature(x="vector"), 
  function(x, nrow, ncol, type="matrix", ..., bldim=.pbd_env$BLDIM, ICTXT=.pbd_env$ICTXT)
  {
    type <- pbdMPI::comm.match.arg(type, c("matrix", "ddmatrix"))
    
    if (missing(nrow) && missing(ncol))
      nrow <- ncol <- length(x)
    else if (missing(nrow) && !missing(ncol))
      nrow <- ncol
    else if (missing(ncol) && !missing(nrow))
      ncol <- nrow
    
    if (type=="ddmatrix")
    {
      if (length(bldim)==1)
        bldim <- rep(bldim, 2)
      
      dim <- c(nrow, ncol)
      ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT)
      
      descx <- base.descinit(dim=dim, bldim=bldim, ldim=ldim, ICTXT=ICTXT)
      
      out <- base.ddiagmk(diag=x, descx=descx)
      ret <- new("ddmatrix", Data=out, dim=dim, ldim=ldim, bldim=bldim, ICTXT=ICTXT)
    }
    else
      ret <- base::diag(x=x, nrow=nrow, ncol=ncol)
    
    return( ret )
  }
)



#' @rdname diag-constructors
#' @export
setMethod("diag", signature(x="character"), 
  function(x, nrow, ncol, type="matrix", ..., min=0, max=1, mean=0, sd=1, rate=1, shape, scale=1, bldim=.pbd_env$BLDIM, ICTXT=.pbd_env$ICTXT)
  {
    type <- pbdMPI::comm.match.arg(type, c("matrix", "ddmatrix"))
    data <- pbdMPI::comm.match.arg(x, c("runif", "uniform", "rnorm", "normal", "rexp", "exponential", "rweibull", "weibull"))
    
    if (missing(nrow) && missing(ncol))
      nrow <- ncol <- length(x)
    else if (missing(nrow) && !missing(ncol))
      nrow <- ncol
    else if (missing(ncol) && !missing(nrow))
      ncol <- nrow
    
    dim <- c(nrow, ncol)
    
    if (type=="ddmatrix")
    {
      if (length(bldim)==1)
        bldim <- rep(bldim, 2)
      
      ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT)
      
      if (!base.ownany(dim=dim, bldim=bldim, ICTXT=ICTXT))
        Data <- matrix(0.0, 1, 1)
      else
      {
        if (data=="runif" || data=="uniform")
          Data <- stats::runif(n=max(ldim), min=min, max=max)
        else if (data=="rnorm" || data=="normal")
          Data <- stats::rnorm(n=max(ldim), mean=mean, sd=sd)
        else if (data=="rexp" || data=="exponential")
          Data <- stats::rexp(n=max(ldim), rate=rate)
        else if (data=="rweibull" || data=="weibull")
          Data <- stats::rweibull(n=max(ldim), shape=shape, scale=scale)
      }
      
      descx <- base.descinit(dim=dim, bldim=bldim, ldim=ldim, ICTXT=ICTXT)
      
      out <- base.ddiagmk(diag=Data, descx=descx)
      ret <- new("ddmatrix", Data=out, dim=dim, ldim=ldim, bldim=bldim, ICTXT=ICTXT)
    }
    else
    {
      if (data=="runif" || data=="uniform")
        Data <- stats::runif(prod(dim), min=min, max=max)
      else if (data=="rnorm" || data=="normal")
        Data <- stats::rnorm(prod(dim), mean=mean, sd=sd)
      else if (data=="rexp" || data=="exponential")
        Data <- stats::rexp(prod(dim), rate=rate)
      else if (data=="rweibull" || data=="weibull")
        Data <- stats::rnorm(prod(dim), min=min, max=max)
      
      ret <- base::diag(x=Data, nrow=nrow, ncol=ncol)
    }
    
    return( ret )
  }
)



#' @rdname diag-constructors
#' @export
setMethod("diag", signature(x="ddmatrix"),
  function(x)
  {
    descx <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
    
    ret <- base.ddiagtk(x=x@Data, descx=descx)
    
    return( ret )
  }
)



# dealing with R being annoying
#' @rdname diag-constructors
#' @export
setMethod("diag", signature(x="matrix"), 
  function(x, nrow, ncol)
    base::diag(x=x)
)
