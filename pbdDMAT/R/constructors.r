#' Distributed Matrix Creation
#' 
#' Methods for simple construction of distributed matrices.
#' 
#' These methods are simplified methods of creating distributed matrices,
#' including random ones.  These methods involve only local computations, i.e.,
#' no communication is performed in the construction of a \code{ddmatrix} using
#' these methods (in contrast to using \code{as.ddmatrix()} et al).
#' 
#' For non-character inputs, the methods attempt to mimic R as closely as
#' possible.  So \code{ddmatrix(1:3, 5, 7)} produces the distributed analogue
#' of \code{matrix(1:3, 5, 7)}.
#' 
#' For character inputs, you may also specify additional parametric family
#' information.
#' 
#' The functions predicated with \code{.local} generate data with a fixed local
#' dimension, i.e., each processor gets an identical amount of data.  Likewise,
#' the remaining functions generate a fixed global amount of data, and each
#' processor may or may not have an identical amount of local data.
#' 
#' To ensure good random number generation, you should only consider using the
#' character methods with the \code{comm.set.seed()} function from pbdMPI which
#' uses the method of L'Ecuyer via the rlecuyer package.
#' 
#' @param data 
#' optional data vector.
#' @param nrow 
#' number of rows.  Global rows for \code{ddmatrix()}. Local rows
#' for \code{ddmatrix.local()}.  See details below.
#' @param ncol 
#' number of columns.  Global columns for \code{ddmatrix()}.  Local
#' columns for \code{ddmatrix.local()}.  See details below.
#' @param byrow 
#' logical. If \code{FALSE} then the distributed matrix will be
#' filled by column major storage, otherwise row-major.
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
#' dx <- ddmatrix(data="rnorm", nrow=5, ncol=6, mean=10, sd=100, bldim=2)
#' dx
#' 
#' dy <- ddmatrix(data=1:4, nrow=7, ncol=5, bldim=2)
#' dy
#' 
#' finalize()
#' }
#' 
#' @keywords Data Generation
#' @name ddmatrix-constructors
#' @rdname ddmatrix-constructors
NULL



#' @rdname ddmatrix-constructors
#' @export
setGeneric(name="ddmatrix", 
  function(data, ...) 
    standardGeneric("ddmatrix"), 
  package="pbdDMAT"
)



### TODO
##' @rdname ddmatrix-constructors
##' @export
#setMethod("dmat", signature(data="missing"), 
#  function(data, nrow=1, ncol=1, byrow=FALSE, ..., bldim=.pbd_env$BLDIM, ICTXT=.pbd_env$ICTXT)
#  {
#    data <- NA
#    ret <- ddmatrix(data=data, nrow=nrow, ncol=ncol, byrow=byrow, bldim=bldim, ICTXT=ICTXT)
#    
#    return( ret )
#  }
#)



#' @rdname ddmatrix-constructors
#' @export
setMethod("ddmatrix", signature(data="ddmatrix"), 
  function(data, nrow=1, ncol=1, byrow=FALSE, ..., bldim=.pbd_env$BLDIM, ICTXT=.pbd_env$ICTXT)
  {
    if (length(bldim)==1)
      bldim <- rep(bldim, 2)
    
    if (nrow==data@dim[1L] && ncol==data@dim[2L])
      return( data )
    else {
      comm.stop("can't do this yet") #FIXME
    }
    
  }
)



#' @rdname ddmatrix-constructors
#' @export
setMethod("ddmatrix", signature(data="missing"), 
  function(data, nrow=1, ncol=1, byrow=FALSE, ..., bldim=.pbd_env$BLDIM, ICTXT=.pbd_env$ICTXT)
  {
    data <- NA
    ret <- ddmatrix(data=data, nrow=nrow, ncol=ncol, byrow=byrow, bldim=bldim, ICTXT=ICTXT)
    
    return( ret )
  }
)



#' @rdname ddmatrix-constructors
#' @export
setMethod("ddmatrix", signature(data="vector"), 
  function(data, nrow=1, ncol=1, byrow=FALSE, ..., bldim=.pbd_env$BLDIM, ICTXT=.pbd_env$ICTXT)
  {
    if (nrow < 1)
      comm.stop("invalid 'nrow'")
    if (ncol < 1)
      comm.stop("invalid 'ncol'")
    
    if (length(bldim)==1)
      bldim <- rep(bldim, 2)
    
    if (missing(nrow))
      nrow <- 1
    if (missing(ncol))
      ncol <- 1
    
    ldata <- base::length(data)
    
    if (ldata > 1){
      if (nrow==1){
        if (ncol==1)
          nrow <- ldata
        else {
          nrow <- ceiling(ldata / ncol)
        }
      }
      else if (ncol==1){
        ncol <- ceiling(ldata / nrow)
      }
      
      dim <- c(nrow, ncol)
      ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT)
      
      Data <- matrix(0.0, ldim[1L], ldim[2L])
      
      descx <- base.descinit(dim=dim, bldim=bldim, ldim=ldim, ICTXT=ICTXT)
      
      MARGIN <- as.integer(byrow) + 1L
      
      Data <- base.pdsweep(x=Data, descx=descx, vec=data, MARGIN=MARGIN, FUN="+")
    } 
    else {
      dim <- c(nrow, ncol)
      ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT)
      
      if (!base.ownany(dim=dim, bldim=bldim, ICTXT=ICTXT))
        Data <- matrix(0.0, 1, 1)
      else
        Data <- matrix(data, ldim[1L], ldim[2L])
    }
    
    # return
    dx <- new("ddmatrix", Data=Data, dim=dim, ldim=ldim, bldim=bldim, ICTXT=ICTXT)
    
    return( dx )
  }
)



#' @rdname ddmatrix-constructors
#' @export
setMethod("ddmatrix", signature(data="matrix"), 
  function(data, nrow=1, ncol=1, byrow=FALSE, ..., bldim=.pbd_env$BLDIM, ICTXT=.pbd_env$ICTXT)
  {
    dim(data) <- NULL
    ret <- ddmatrix(data=data, nrow=nrow, ncol=ncol, byrow=byrow, bldim=bldim, ICTXT=ICTXT)
    
    return( ret )
  }
)



#' @rdname ddmatrix-constructors
#' @export
setMethod("ddmatrix", signature(data="character"), 
  function(data, nrow=1, ncol=1, byrow=FALSE, ..., min=0, max=1, mean=0, sd=1, rate=1, shape, scale=1, bldim=.pbd_env$BLDIM, ICTXT=.pbd_env$ICTXT)
  {
    data <- pbdMPI::comm.match.arg(data, c("runif", "uniform", "rnorm", "normal", "rexp", "exponential", "rweibull", "weibull"))
    
    if (length(bldim)==1)
      bldim <- rep(bldim, 2)
    
    dim <- c(nrow, ncol)
    ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT)
    
    if (!base.ownany(dim=dim, bldim=bldim, ICTXT=ICTXT))
      Data <- matrix(0.0, 1, 1)
    else {
      if (data=="runif" || data=="uniform")
        Data <- stats::runif(n=prod(ldim), min=min, max=max)
      else if (data=="rnorm" || data=="normal")
        Data <- stats::rnorm(n=prod(ldim), mean=mean, sd=sd)
      else if (data=="rexp" || data=="exponential")
        Data <- stats::rexp(n=prod(ldim), rate=rate)
      else if (data=="rweibull" || data=="weibull")
        Data <- stats::rweibull(n=prod(ldim), shape=shape, scale=scale)
      
      dim(Data) <- ldim
    }
    
    dx <- new("ddmatrix", Data=Data, dim=dim, ldim=ldim, bldim=bldim, ICTXT=ICTXT)
    
    return( dx )
  }
)



### local versions; not sure how useful this is to anyone, but why not?
#' @rdname ddmatrix-constructors
#' @export
setGeneric(name="ddmatrix.local", 
  function(data, ...) 
    standardGeneric("ddmatrix.local"), 
  package="pbdDMAT"
)



#' @rdname ddmatrix-constructors
#' @export
setMethod("ddmatrix.local", signature(data="missing"), 
  function(data, nrow=1, ncol=1, byrow=FALSE, ..., bldim=.pbd_env$BLDIM, ICTXT=.pbd_env$BLDIM)
  {
    data <- NA
    ret <- ddmatrix.local(data=data, nrow=nrow, ncol=ncol, byrow=byrow, bldim=bldim, ICTXT=ICTXT)
    
    return( ret )
  }
)



#' @rdname ddmatrix-constructors
#' @export
setMethod("ddmatrix.local", signature(data="vector"), 
  function(data, nrow=1, ncol=1, byrow=FALSE, ..., bldim=.pbd_env$BLDIM, ICTXT=.pbd_env$ICTXT)
  {
    if (nrow < 1)
      comm.stop("invalid 'nrow'")
    if (ncol < 1)
      comm.stop("invalid 'ncol'")
    
    if (length(bldim)==1)
      bldim <- rep(bldim, 2)
    
    ldim <- c(nrow, ncol)
    
    blacs_ <- base.blacs(ICTXT=ICTXT)
    nprows <- blacs_$NPROW
    npcols <- blacs_$NPCOL
    
    dim <- c(nprows*ldim[1L], npcols*ldim[2L])
    
    # bldim
    if (any( (dim %% bldim) != 0 )){
      comm.warning("at least one margin of 'bldim' does not divide the global dimension.\n")
      
      bldim[1L] <- base.nbd(ldim[1L], bldim[1L])
      bldim[2L] <- base.nbd(ldim[2L], bldim[2L])
      comm.cat(paste("Using bldim of ", bldim[1L], "x", bldim[2L], "\n\n", sep=""), quiet=TRUE)
    }
    
    if (length(data) > 1){
      Data <- matrix(0.0, ldim[1L], ldim[2L])
      
      descx <- base.descinit(dim=dim, bldim=bldim, ldim=ldim, ICTXT=ICTXT)
      
      MARGIN <- as.integer(byrow) + 1L
      
      Data <- base.pdsweep(x=Data, descx=descx, vec=data, MARGIN=MARGIN, FUN="+")
    } 
    else {
      if (!base.ownany(dim=dim, bldim=bldim, ICTXT=ICTXT))
        Data <- matrix(0.0, 1, 1)
      else
        Data <- matrix(data, ldim[1L], ldim[2L])
    }
    
    # return
    dx <- new("ddmatrix", Data=Data, dim=dim, ldim=ldim, bldim=bldim, ICTXT=ICTXT)
    
    return( dx )
  }
)



#' @rdname ddmatrix-constructors
#' @export
setMethod("ddmatrix.local", signature(data="matrix"), 
  function(data, nrow=1, ncol=1, byrow=FALSE, ..., bldim=.pbd_env$BLDIM, ICTXT=.pbd_env$ICTXT)
  {
    dim(data) <- NULL
    ret <- ddmatrix.local(data=data, nrow=nrow, ncol=ncol, byrow=byrow, bldim=bldim, ICTXT=ICTXT)
    
    return( ret )
  }
)



#' @rdname ddmatrix-constructors
#' @export
setMethod("ddmatrix.local", signature(data="character"), 
  function(data, nrow=1, ncol=1, byrow=FALSE, ..., min=0, max=1, mean=0, sd=1, rate=1, shape, scale=1, bldim=.pbd_env$BLDIM, ICTXT=.pbd_env$ICTXT)
  {
    if (nrow < 1)
      comm.stop("invalid 'nrow'")
    if (ncol < 1)
      comm.stop("invalid 'ncol'")
    
    if (length(bldim)==1)
      bldim <- rep(bldim, 2)
    
    data <- pbdMPI::comm.match.arg(data, c("runif", "uniform", "rnorm", "normal", "rexp", "exponential", "rweibull", "weibull"))
    
    ldim <- c(nrow, ncol)
    
    blacs_ <- base.blacs(ICTXT=ICTXT)
    nprows <- blacs_$NPROW
    npcols <- blacs_$NPCOL
    
    dim <- c(nprows*ldim[1L], npcols*ldim[2L])
    
    # bldim
    if (any( (dim %% bldim) != 0 )){
      comm.warning("at least one margin of 'bldim' does not divide the global dimension.\n")
      
      bldim[1L] <- base.nbd(ldim[1L], bldim[1L])
      bldim[2L] <- base.nbd(ldim[2L], bldim[2L])
      comm.cat(paste("Using bldim of ", bldim[1L], "x", bldim[2L], "\n\n", sep=""), quiet=TRUE)
    }
    
    
    if (!base.ownany(dim=dim, bldim=bldim, ICTXT=ICTXT))
      Data <- matrix(0.0, 1, 1)
    else {
      if (data=="runif" || data=="uniform")
        Data <- matrix(runif(prod(ldim), min=min, max=max), ldim[1L], ldim[2L])
      else if (data=="rnorm" || data=="normal")
        Data <- matrix(rnorm(prod(ldim), mean=mean, sd=sd), ldim[1L], ldim[2L])
      else if (data=="rexp" || data=="exponential")
        Data <- matrix(rexp(prod(ldim), rate=rate), ldim[1L], ldim[2L])
      else if (data=="rweibull" || data=="weibull")
        Data <- matrix(rweibull(n=prod(ldim), shape=shape, scale=scale), ldim[1L], ldim[2L])
    }
    
    dx <- new("ddmatrix", Data=Data, dim=dim, ldim=ldim, bldim=bldim, ICTXT=ICTXT)
    
    return( dx )
  }
)
