#' QR Decomposition Methods
#' 
#' \code{qr()} takes the QR decomposition.
#' 
#' \code{qr.Q()} recovers Q from the output of \code{qr()}.
#' 
#' \code{qr.R()} recovers R from the output of \code{qr()}.
#' 
#' \code{qr.qy()} multiplies \code{y} by Q.
#' 
#' \code{qr.qty()} multiplies \code{y} by the transpose of Q.
#' 
#' Functions for forming a QR decomposition and for using the outputs of these
#' numerical QR routines.
#' 
#' @param x,y 
#' numeric distributed matrices for \code{qr()}. Otherwise, \code{x}
#' is a list, namely the return from \code{qr()}.
#' @param tol 
#' logical value, determines whether or not columns are zero
#' centered.
#' @param complete 
#' logical expression of length 1.  Indicates whether an
#' arbitrary orthogonal completion of the Q or X matrices is to be made, or
#' whether the R matrix is to be completed by binding zero-value rows beneath
#' the square upper triangle.
#' @param Dvec 
#' Not implemented for objects of class \code{ddmatrix}.  vector
#' (not matrix) of diagonal values.  Each column of the returned Q will be
#' multiplied by the corresponding diagonal value.  Defaults to all 1's.
#' 
#' @return 
#' \code{qr()} returns a list consisting of: \code{qr} - \code{rank} -
#' calculated numerical rank, \code{tau} - \code{pivot} - \code{"class"} -
#' attribute "qr".
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
#' # don't do this in production code
#' x <- matrix(1:9, 3)
#' x <- as.ddmatrix(x)
#' 
#' Q <- qr.Q(qr(x))
#' print(Q)
#' 
#' finalize()
#' }
#' 
#' @keywords Methods Linear Algebra
#' @aliases qr.Q qr.R qr.qty qr.qy
#' @name qr
#' @rdname qr
NULL



setGeneric(name="qr", 
  function(x, ...)
    standardGeneric("qr"),
  package="pbdDMAT"
)

setGeneric(name="qr.Q", 
  function(x, ...)
    standardGeneric("qr.Q"),
  package="pbdDMAT"
)

setGeneric(name="qr.R", 
  function(x, ...)
    standardGeneric("qr.R"),
  package="pbdDMAT"
)

setGeneric(name="qr.qy", 
  function(x, ...)
    standardGeneric("qr.qy"),
  package="pbdDMAT"
)

setGeneric(name="qr.qty", 
  function(x, ...)
    standardGeneric("qr.qty"),
  package="pbdDMAT"
)



#' @rdname qr
#' @export
setMethod("qr", signature(x="ddmatrix"), 
  function (x, tol = 1e-07)
  {
    # Matrix descriptors
    descx <- base.descinit(x@dim, x@bldim, x@ldim, ICTXT=x@ICTXT)
    
    m <- descx[3L]
    n <- descx[4L]
    
    # qr
    out <- base.rpdgeqpf(tol=tol, m=m, n=n, x=x@Data, descx=descx)
    
    # make sure processors who own nothing have a real value (not a null
    # pointer) for the numerical rank
#    if (comm.rank()!=0)
#      rank <- 0
#    else
#      rank <- out$rank
#    
#    rank <- allreduce(rank)
    
    
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT)){
      x@Data <- out$qr
    }
    
    ret <- list(qr=x, rank=out$rank, tau=out$tau, pivot=out$pivot)
    
    attr(ret, "class") <- "qr"
    
    return( ret )
  }
)



#' @rdname qr
#' @export
setMethod("qr.Q", signature(x="ANY"), 
  function(x, complete = FALSE,  Dvec)
    {
      # x is of class qr
      
      if (is.ddmatrix(x$qr)){
        # complete/Dvec options
        qr <- x$qr
        
        if (qr@dim[1L] < qr@dim[2L])
          qr <- qr[, 1L:x$rank]
        
        # Matrix descriptors
        descqr <- base.descinit(qr@dim, qr@bldim, qr@ldim, ICTXT=qr@ICTXT)
        
        m <- descqr[3]
        n <- descqr[4]
        
        k <- x$rank
        
        ret <- base.rpdorgqr(m=m, n=n, k=k, qr=qr@Data, descqr=descqr, tau=x$tau)
        
        if (base.ownany(dim=qr@dim, bldim=qr@bldim, ICTXT=qr@ICTXT)){
          qr@Data <- ret
        }
        
        return( qr )
        
      } else {
        dqr <- dim(x$qr)
        cmplx <- mode(x$qr) == "complex"
        ret <- base::qr.Q(qr=x, complete=complete, Dvec=Dvec)
      }
      
      return( ret )
  }
)



# qr.R
dmat.qr.R <- function(qr, complete=FALSE)
{
  ret <- qr$qr
  
  if (!complete){
    if (min(ret@dim)!=ret@dim[1L])
      ret <- ret[1L:min(ret@dim), ]
  }
  
  descx <- base.descinit(dim=ret@dim, bldim=ret@bldim, ldim=ret@ldim, ICTXT=ret@ICTXT)
  
  ret@Data <- base.tri2zero(x=ret@Data, descx=descx, uplo='L', diag='N')
  
  # not particularly efficient, but I don't expect this to get any real use...
  rank <- qr$rank
  n <- ret@dim[1L]
  p <- ret@dim[2L]
  mn <- min(ret@dim)
  
  if (rank < p){
    if (n>p)
      for (i in (rank+1L):mn)
        ret[i,i] <- 0
  }
  
  return(ret)
}



#' @rdname qr
#' @export
setMethod("qr.R", signature(x="ANY"), 
  function(x, complete = FALSE) 
  {
    qr <- x
    
    if (is.ddmatrix(qr$qr)){
      ret <- dmat.qr.R(qr=qr, complete=complete)
    } else {
      ret <- base::qr.R(qr=qr, complete=complete)
    }
    
    return( ret )
  }
)



#' @rdname qr
#' @export
setMethod("qr.qy", signature(x="ANY"), 
  function(x, y)
  {
    if (is.ddmatrix(x$qr)){
      
      qr <- x$qr
      
      # Matrix descriptors
      descqr <- base.descinit(qr@dim, qr@bldim, qr@ldim, ICTXT=qr@ICTXT)
      descc <- base.descinit(y@dim, y@bldim, y@ldim, ICTXT=y@ICTXT)
      
      m <- descqr[3L]
      n <- y@dim[2L]
      k <- x$rank
      
      out <- base.rpdormqr(side='L', trans='N', m=m, n=n, k=k, qr=qr@Data, descqr=descqr, tau=x$tau, c=y@Data, descc=descc)
      
      y@Data <- out
      
      return(y)
      
    } else {
      ret <- base::qr.qy(qr=x, y=y)
    }
    
    return( ret )
  }
)



#' @rdname qr
#' @export
setMethod("qr.qty", signature(x="ANY"), 
  function(x, y)
  {
    if (is.ddmatrix(x$qr)){
      
      qr <- x$qr
      
      # Matrix descriptors
      descqr <- base.descinit(qr@dim, qr@bldim, qr@ldim, ICTXT=qr@ICTXT)
      descc <- base.descinit(y@dim, y@bldim, y@ldim, ICTXT=y@ICTXT)
      
      m <- descqr[3L]
      n <- y@dim[2L]
      k <- x$rank
      
      out <- base.rpdormqr(side='L', trans='T', m=m, n=n, k=k, qr=qr@Data, descqr=descqr, tau=x$tau, c=y@Data, descc=descc)
      
      y@Data <- out
      
      return(y)
      
    } else {
      ret <- base::qr.qty(qr=x, y=y)
    }
    
    return( ret )
  }
)





