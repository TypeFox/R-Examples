#' Logical Comparisons
#' 
#' Logical comparisons.
#' 
#' Performs the indicated logical comparison.
#' 
#' If \code{na.rm} is \code{TRUE} and only \code{NA}'s are present, then
#' \code{TRUE} is returned.
#' 
#' @param e1,e2,x
#' distributed matrix or numeric vector
#' @param na.rm 
#' logical, indicating whether or not \code{NA}'s should first be
#' removed. If not and an NA is present, \code{NA} is returned.
#' 
#' @return 
#' Returns a distributed matrix.
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
#' x <- matrix(sample(0, 1, 9, replace=T), 3)
#' comm.print(x)
#' 
#' x <- as.ddmatrix(x, bldim=2)
#' 
#' y <- any(x)
#' comm.print(y)
#' 
#' finalize()
#' }
#' 
#' @keywords Methods Extraction Type
#' @name Comparators
#' @rdname Comparators
NULL


# -------------------
# ddmatrix Comparators
# -------------------

#' @rdname Comparators
#' @export
setMethod("==", signature(e1="ddmatrix", e2="ddmatrix"),
  function(e1, e2)
  {
    base.checkem(x=e1, y=e2, checks=1:3)
    if (base.ownany(dim=e1@dim, bldim=e1@bldim, ICTXT=e1@ICTXT))
      e1@Data <- e1@Data == e2@Data
    
    return( e1 )
  }
) 

#' @rdname Comparators
#' @export
setMethod("!=", signature(e1="ddmatrix", e2="ddmatrix"),
  function(e1, e2)
  {
    base.checkem(x=e1, y=e2, checks=1:3)
    if (base.ownany(dim=e1@dim, bldim=e1@bldim, ICTXT=e1@ICTXT))
      e1@Data <- e1@Data != e2@Data
    
    return( e1 )
  }
) 

#' @rdname Comparators
#' @export
setMethod("all", signature(x="ddmatrix"),
  function(x, na.rm=FALSE)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      ret <- base::all(x@Data)
    else
      ret <- 1
    
    ret <- as.logical( pbdMPI::allreduce(ret, op='min') )
    
    return(ret)
  }
) 

#' @rdname Comparators
#' @export
setMethod("any", signature(x="ddmatrix"),
  function(x, na.rm=FALSE)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      ret <- base::any(x@Data)
    else
      ret <- 0
    
    ret <- as.logical( pbdMPI::allreduce(ret, op='max') )
    
    return(ret)
  }
) 

#' @rdname Comparators
#' @export
setMethod("<", signature(e1="ddmatrix", e2="ddmatrix"),
  function(e1, e2)
  {
    base.checkem(x=e1, y=e2, checks=1:3)
    if (base.ownany(dim=e1@dim, bldim=e1@bldim, ICTXT=e1@ICTXT))
      e1@Data <- e1@Data < e2@Data
    
    return(e1)
  }
) 

#' @rdname Comparators
#' @export
setMethod(">", signature(e1="ddmatrix", e2="ddmatrix"),
  function(e1, e2)
  {
    base.checkem(x=e1, y=e2, checks=1:3)
    if (base.ownany(dim=e1@dim, bldim=e1@bldim, ICTXT=e1@ICTXT))
      e1@Data <- e1@Data > e2@Data
    
    return(e1)
  }
) 

#' @rdname Comparators
#' @export
setMethod("<=", signature(e1="ddmatrix", e2="ddmatrix"),
  function(e1, e2)
  {
    base.checkem(x=e1, y=e2, checks=1:3)
    if (base.ownany(dim=e1@dim, bldim=e1@bldim, ICTXT=e1@ICTXT))
      e1@Data <- e1@Data <= e2@Data
    
    return(e1)
  }
) 

#' @rdname Comparators
#' @export
setMethod(">=", signature(e1="ddmatrix", e2="ddmatrix"),
  function(e1, e2)
  {
    base.checkem(x=e1, y=e2, checks=1:3)
    if (base.ownany(dim=e1@dim, bldim=e1@bldim, ICTXT=e1@ICTXT))
      e1@Data <- e1@Data >= e2@Data
    
    return(e1)
  }
) 

#' @rdname Comparators
#' @export
setMethod("|", signature(e1="ddmatrix", e2="ddmatrix"),
  function(e1, e2)
  {
    base.checkem(x=e1, y=e2, checks=1:3)
    if (base.ownany(dim=e1@dim, bldim=e1@bldim, ICTXT=e1@ICTXT))
      e1@Data <- e1@Data | e2@Data
    
    return( e1 )
  }
) 

#' @rdname Comparators
#' @export
setMethod("&", signature(e1="ddmatrix", e2="ddmatrix"),
  function(e1, e2)
  {
    base.checkem(x=e1, y=e2, checks=1:3)
    if (base.ownany(dim=e1@dim, bldim=e1@bldim, ICTXT=e1@ICTXT))
      e1@Data <- e1@Data & e2@Data
    
    return( e1 )
  }
) 


# -------------------
# ddmatrix-vector Comparators
# -------------------

#' @rdname Comparators
#' @export
setMethod("<", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- base::length(e2)
    if ( (prod(dim)%%len > 0 && len%%prod(dim) > 0) && len > 1)
      comm.warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, ICTXT=e1@ICTXT)){
      if (len==1)
        e1@Data <- e1@Data<e2
      else {
        descx <- base.descinit(dim=e1@dim, bldim=e1@bldim, ldim=e1@ldim, ICTXT=e1@ICTXT)
        out <- base.rl2blas(x=e1@Data, descx=descx, vec=e2, FUN=7)
        
        dim(out) <- e1@ldim
        if (!is.logical(out))
          storage.mode(out) <- "logical"
        
        e1@Data <- out
      }
    }
    
    return(e1)
  }
)

#' @rdname Comparators
#' @export
setMethod("<", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2)
    e2>e1
)

#' @rdname Comparators
#' @export
setMethod(">", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- base::length(e2)
    if ( (prod(dim)%%len > 0 && len%%prod(dim) > 0) && len > 1)
      comm.warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, ICTXT=e1@ICTXT)){
      if (len==1)
        e1@Data <- e1@Data>e2
      else {
        descx <- base.descinit(dim=e1@dim, bldim=e1@bldim, ldim=e1@ldim, ICTXT=e1@ICTXT)
        out <- base.rl2blas(x=e1@Data, descx=descx, vec=e2, FUN=8)
        
        dim(out) <- e1@ldim
        if (!is.logical(out))
          storage.mode(out) <- "logical"
        
        e1@Data <- out
      }
    }
    return(e1)
  }
)

#' @rdname Comparators
#' @export
setMethod(">", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2)
    e2<e1
)

#' @rdname Comparators
#' @export
setMethod("<=", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- base::length(e2)
    if ( (prod(dim)%%len > 0 && len%%prod(dim) > 0) && len > 1)
      comm.warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, ICTXT=e1@ICTXT)){
      if (len==1)
        e1@Data <- e1@Data<=e2
      else {
        descx <- base.descinit(dim=e1@dim, bldim=e1@bldim, ldim=e1@ldim, ICTXT=e1@ICTXT)
        out <- base.rl2blas(x=e1@Data, descx=descx, vec=e2, FUN=9)
        
        dim(out) <- e1@ldim
        if (!is.logical(out))
          storage.mode(out) <- "logical"
        
        e1@Data <- out
      }
    }
    else {
      e1@Data <- matrix(0.0, 1, 1)
    }
    return(e1)
  }
)

#' @rdname Comparators
#' @export
setMethod("<=", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2)
    e2>=e1
)

#' @rdname Comparators
#' @export
setMethod(">=", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- base::length(e2)
    if ( (prod(dim)%%len > 0 && len%%prod(dim) > 0) && len > 1)
      comm.warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, ICTXT=e1@ICTXT)){
      if (len==1)
        e1@Data <- e1@Data>=e2
      else {
        descx <- base.descinit(dim=e1@dim, bldim=e1@bldim, ldim=e1@ldim, ICTXT=e1@ICTXT)
        out <- base.rl2blas(x=e1@Data, descx=descx, vec=e2, FUN=10)
        
        dim(out) <- e1@ldim
        if (!is.logical(out))
          storage.mode(out) <- "logical"
        
        e1@Data <- out
      }
    }
    return(e1)
  }
)

#' @rdname Comparators
#' @export
setMethod(">=", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2)
    e2<=e1
)

#' @rdname Comparators
#' @export
setMethod("==", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- base::length(e2)
    if ( (prod(dim)%%len > 0 && len%%prod(dim) > 0) && len > 1)
      comm.warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, ICTXT=e1@ICTXT)){
      if (len==1)
        e1@Data <- e1@Data==e2
      else {
        descx <- base.descinit(dim=e1@dim, bldim=e1@bldim, ldim=e1@ldim, ICTXT=e1@ICTXT)
        out <- base.rl2blas(x=e1@Data, descx=descx, vec=e2, FUN=11)
        
        dim(out) <- e1@ldim
        if (!is.logical(out))
          storage.mode(out) <- "logical"
        
        e1@Data <- out
      }
    }
    return(e1)
  }
)

#' @rdname Comparators
#' @export
setMethod("==", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2)
    e2==e1
)

#' @rdname Comparators
#' @export
setMethod("!=", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- base::length(e2)
    if ( (prod(dim)%%len > 0 && len%%prod(dim) > 0) && len > 1)
      comm.warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, ICTXT=e1@ICTXT)){
      if (len==1)
        e1@Data <- e1@Data != e2
      else {
        descx <- base.descinit(dim=e1@dim, bldim=e1@bldim, ldim=e1@ldim, ICTXT=e1@ICTXT)
        out <- base.rl2blas(x=e1@Data, descx=descx, vec=e2, FUN=11)
        
        dim(out) <- e1@ldim
        if (!is.logical(out))
          storage.mode(out) <- "logical"
        
        out <- !out
        
        e1@Data <- out
      }
    }
    return(e1)
  }
)

#' @rdname Comparators
#' @export
setMethod("!=", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2)
    e2!=e1
)

#' @rdname Comparators
#' @export
setMethod("|", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- base::length(e2)
    if ( (prod(dim)%%len > 0 && len%%prod(dim) > 0) && len > 1)
      comm.warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, ICTXT=e1@ICTXT)){
      if (len==1)
        e1@Data <- e1@Data | e2
      else {
        descx <- base.descinit(dim=e1@dim, bldim=e1@bldim, ldim=e1@ldim, ICTXT=e1@ICTXT)
        out <- base.rl2blas(x=e1@Data, descx=descx, vec=e2, FUN=0)
        
        dim(out) <- e1@ldim
        if (!is.logical(out))
          storage.mode(out) <- "logical"
        
        e1@Data <- out
      }
    }
    return(e1)
  }
)

#' @rdname Comparators
#' @export
setMethod("|", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2)
    e2 | e1
)

#' @rdname Comparators
#' @export
setMethod("&", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- base::length(e2)
    if ( (prod(dim)%%len > 0 && len%%prod(dim) > 0) && len > 1)
      comm.warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, ICTXT=e1@ICTXT)){
      if (len==1)
        e1@Data <- e1@Data & e2
      else {
        descx <- base.descinit(dim=e1@dim, bldim=e1@bldim, ldim=e1@ldim, ICTXT=e1@ICTXT)
        out <- base.rl2blas(x=e1@Data, descx=descx, vec=e2, FUN=2)
        
        dim(out) <- e1@ldim
        if (!is.logical(out))
          storage.mode(out) <- "logical"
        
        e1@Data <- out
      }
    }
    return(e1)
  }
)

#' @rdname Comparators
#' @export
setMethod("&", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2)
    e2 & e1
)


