#' Sweep
#' 
#' Sweep vector or ddmatrix from a distributed matrix.
#' 
#' @param x 
#' numeric distributed matrix.
#' @param MARGIN 
#' subscript over which the function will be applied
#' @param STATS 
#' array to be swept out.
#' @param FUN 
#' function used in the sweep. Only \code{+}, \code{-}, \code{*},
#' and \code{/} are accepted.  For more general operations, use \code{apply()}.
#' @param check.margin,...
#' Ignored.
#' 
#' @return 
#' Returns a distributed matrix.
#' 
#' @keywords Methods
#' 
#' @name sweep
#' @rdname sweep
NULL



#' @rdname sweep
#' @export
setGeneric(name = "sweep", useAsDefault = base::sweep, package="pbdDMAT")

#' @rdname sweep
#' @export
setMethod("sweep", signature(x="ddmatrix", STATS="vector"),
  function(x, MARGIN, STATS, FUN = "-")
  {
    # checks
    if ( !(FUN %in% c("+", "-", "*", "/")) )
      comm.stop("Error : invalid argument 'FUN'")
    
    if (MARGIN != 1 && MARGIN != 2)
      comm.stop("Error : argument 'MARGIN' must be 1 or 2")
    
    if ( is.matrix(STATS) )
      dim(STATS) <- NULL
    
    descx <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
    
    ret <- base.pdsweep(x=x@Data, descx=descx, vec=STATS, MARGIN=MARGIN, FUN=FUN)
    
    x@Data <- ret
    
    return( x )
  }
)

# trapped in bad code factory, send help
#' @rdname sweep
#' @export
setMethod("sweep", signature(x="ddmatrix", STATS="ddmatrix"),
  function(x, MARGIN, STATS, FUN = "-")
  {
    # checks
    if ( !(FUN %in% c("+", "-", "*", "/")) )
      comm.stop("Error : invalid argument 'FUN'")
    
    if (MARGIN != 1L && MARGIN != 2L)
      comm.stop("Error : argument 'MARGIN' must be 1 or 2")
    
    if (any(x@bldim != STATS@bldim))
      comm.stop("Error : blocking dimensions of 'x' and 'STATS' must be identical")
    
    if (x@ICTXT != STATS@ICTXT)
      comm.stop("Error : ICTXT of 'x' and 'STATS' must be the same")
    
    # work in place if possible, otherwise cast as global vector to preserve R-like BLAS
    if (all(x@dim==STATS@dim)){
      if (MARGIN == 1)
        ret <- x - STATS
      else if (any(x@dim) == 1)
        ret <- x - STATS
      else {
        # FIXME
        STATS <- as.vector(STATS)
        return( sweep(x=x, MARGIN=MARGIN, STATS=STATS, FUN=FUN) )
      }
    }
    else if (x@dim[MARGIN] != 1){
      if (STATS@dim[2L/MARGIN]==1 && STATS@dim[MARGIN]==x@dim[MARGIN]){
#        x@Data <- base::sweep(x=x@Data, STATS=STATS@Data, MARGIN=MARGIN, FUN=FUN)
        vec <- STATS@Data
        
        len <- dmat.allrowreduce(x=base::length(vec), op='max', ICTXT=x@ICTXT)
        
        if (!base.ownany(dim=STATS@dim, bldim=STATS@bldim, ICTXT=STATS@ICTXT))
          vec <- numeric(len)
        
        vec <- dmat.allrowreduce(x=vec, op='sum', ICTXT=x@ICTXT)
        vec <- vec + 0.0
        
        out <- base::sweep(x=x@Data, STATS=vec, MARGIN=MARGIN, FUN=FUN)
        
        x@Data <- out
        return( x )
      }
      else if (STATS@dim[MARGIN]==1 && STATS@dim[2L/MARGIN]==x@dim[MARGIN]) {
        STATS <- t(STATS)
        vec <- STATS@Data
        dim(vec) <- NULL
        
        len <- dmat.allcolreduce(x=base::length(vec), op='max', ICTXT=x@ICTXT)
        
        if (!base.ownany(dim=STATS@dim, bldim=STATS@bldim, ICTXT=STATS@ICTXT))
          vec <- numeric(len)
        
        vec <- dmat.allcolreduce(x=vec, op='sum', ICTXT=x@ICTXT)
        vec <- vec + 0.0
        
        out <- base::sweep(x=x@Data, STATS=vec, MARGIN=MARGIN, FUN=FUN)
        
        x@Data <- out
        return( x )
      }
      else {
        # FIXME
        comm.print("cast as vec")
        STATS <- as.vector(STATS)
        return( sweep(x=x, MARGIN=MARGIN, STATS=STATS, FUN=FUN) )
      }
    }
    else {
      # FIXME
      STATS <- as.vector(STATS)
      return( sweep(x=x, MARGIN=MARGIN, STATS=STATS, FUN=FUN) )
    }
    
    return( ret )
  }
)

