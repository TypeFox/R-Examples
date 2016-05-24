#' Arithmetic Operators
#' 
#' Binary operations for distributed matrix/distributed matrix and distributed
#' matrix/vector operations.
#' 
#' If \code{e1} and \code{e2} are distributed matrices, then they must be
#' conformable, on the same BLACS context, and have the same blocking
#' dimension.
#' 
#' @param e1,e2
#' numeric distributed matrices or numeric vectors
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
#' x <- ddmatrix(1:9, 3, bldim=2)
#' 
#' y <- (2*x) - x^(.5)
#' print(y)
#' 
#' finalize()
#' }
#' 
#' @name arithmetic
#' @rdname arithmetic
#' @keywords Methods
NULL


# ----------------
# +
# ----------------

#' @rdname arithmetic
#' @export
setMethod("+", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- length(e2)
    if ( (prod(dim)%%len > 0 && len%%prod(dim) > 0) && len > 1)
      comm.warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, ICTXT=e1@ICTXT)){
      if (len==1)
        e1@Data <- e1@Data+e2
      else {
        descx <- base.descinit(dim=e1@dim, bldim=e1@bldim, ldim=e1@ldim, ICTXT=e1@ICTXT)
        e1@Data <- base.rl2blas(x=e1@Data, descx=descx, vec=e2, FUN=0) # FUN=0 for "+"
      }
    }
    
    return(e1)
  }
)



#' @rdname arithmetic
#' @export
setMethod("+", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2)
    e2+e1
)



#' @rdname arithmetic
#' @export
setMethod("+", signature(e1="ddmatrix", e2="ddmatrix"), 
  function(e1, e2){
    base.checkem(x=e1, y=e2, checks=1:3)
    if (base.ownany(dim=e1@dim, bldim=e1@bldim, ICTXT=e1@ICTXT))
      e1@Data <- e1@Data + e2@Data
    
    return(e1)
  }
)



# ----------------
# -
# ----------------

#' @rdname arithmetic
#' @export
setMethod("-", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- length(e2)
    if ( (dim[1]%%len > 0 && len%%dim[1] > 0) && len > 1 )
      comm.warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, ICTXT=e1@ICTXT))
      if (len==1)
        e1@Data <- e1@Data-e2
      else {
        descx <- base.descinit(dim=e1@dim, bldim=e1@bldim, ldim=e1@ldim, ICTXT=e1@ICTXT)
        e1@Data <- base.rl2blas(x=e1@Data, descx=descx, vec=e2, FUN=1) # FUN=1 for "-"
      }
    
    return(e1)
  }
)



#' @rdname arithmetic
#' @export
setMethod("-", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2){
    e2@Data <- -e2@Data
    
    return(e2+e1)
  }
)



#' @rdname arithmetic
#' @export
setMethod("-", signature(e1="ddmatrix", e2="ddmatrix"), 
  function(e1, e2){
    base.checkem(x=e1, y=e2, checks=1:3)
    if (base.ownany(dim=e1@dim, bldim=e1@bldim, ICTXT=e1@ICTXT))
      e1@Data <- e1@Data - e2@Data
    
    return(e1)
  }
)



#' @rdname arithmetic
#' @export
setMethod("-", signature(e1="ddmatrix", e2="missing"), 
  function(e1, e2){
    e1@Data <- -e1@Data
    
    return(e1)
  }
)



# ----------------
# *
# ----------------

#' @rdname arithmetic
#' @export
setMethod("*", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- length(e2)
    if ( (dim[1]%%len > 0 && len%%dim[1] > 0) && len > 1)
      comm.warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, ICTXT=e1@ICTXT))
      if (len==1)
        e1@Data <- e1@Data*e2
      else {
        descx <- base.descinit(dim=e1@dim, bldim=e1@bldim, ldim=e1@ldim, ICTXT=e1@ICTXT)
        e1@Data <- base.rl2blas(x=e1@Data, descx=descx, vec=e2, FUN=2) # FUN=2 for "*"
      }
    
    return(e1)
  }
)



#' @rdname arithmetic
#' @export
setMethod("*", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2)
    return(e2*e1)
)



#' @rdname arithmetic
#' @export
setMethod("*", signature(e1="ddmatrix", e2="ddmatrix"), 
  function(e1, e2){
    base.checkem(x=e1, y=e2, checks=1:3)
    if (base.ownany(dim=e1@dim, bldim=e1@bldim, ICTXT=e1@ICTXT))
      e1@Data <- e1@Data * e2@Data
    
    return(e1)
  }
)



# ----------------
# /
# ----------------

#' @rdname arithmetic
#' @export
setMethod("/", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- length(e2)
    if ( (dim[1]%%len > 0 && len%%dim[1] > 0) && len > 1)
      comm.warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, ICTXT=e1@ICTXT))
      if (len==1)
        e1@Data <- e1@Data/e2
      else {
        descx <- base.descinit(dim=e1@dim, bldim=e1@bldim, ldim=e1@ldim, ICTXT=e1@ICTXT)
        e1@Data <- base.rl2blas(x=e1@Data, descx=descx, vec=e2, FUN=3) # FUN=3 for "/"
      }
      
    return(e1)
  }
)



#' @rdname arithmetic
#' @export
setMethod("/", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2){
    if (base.ownany(dim=e2@dim, bldim=e2@bldim, ICTXT=e2@ICTXT))
      e2@Data <- 1 / e2@Data
    
    return(e2*e1)
  }
)



#' @rdname arithmetic
#' @export
setMethod("/", signature(e1="ddmatrix", e2="ddmatrix"), 
  function(e1, e2){
    base.checkem(x=e1, y=e2, checks=1:3)
    if (base.ownany(dim=e1@dim, bldim=e1@bldim, ICTXT=e1@ICTXT))
      e1@Data <- e1@Data / e2@Data
    
    return(e1)
  }
)



# ----------------
# ^
# ----------------

#' @rdname arithmetic
#' @export
setMethod("^", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- length(e2)
    if ( (dim[1]%%len > 0 && len%%dim[1] > 0) && len > 1)
      comm.warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, ICTXT=e1@ICTXT))
      if (len==1)
        e1@Data <- e1@Data^e2
      else {
        descx <- base.descinit(dim=e1@dim, bldim=e1@bldim, ldim=e1@ldim, ICTXT=e1@ICTXT)
        e1@Data <- base.rl2blas(x=e1@Data, descx=descx, vec=e2, FUN=4) # FUN=4 for "^"
      }
    
    return(e1)
  }
)



## numeric ^ ddmatrix
#setMethod("^", signature(e1="numeric", e2="ddmatrix"), 
#  function(e1, e2){
#    
#  }
#)



#' @rdname arithmetic
#' @export
setMethod("^", signature(e1="ddmatrix", e2="ddmatrix"), 
  function(e1, e2){
    base.checkem(x=e1, y=e2, checks=1:3)
    if (base.ownany(dim=e1@dim, bldim=e1@bldim, ICTXT=e1@ICTXT))
      e1@Data <- e1@Data ^ e2@Data
    
    return(e1)
  }
)



# ----------------
# Modulo stuff --- pretty useless, really
# ----------------

#' @rdname arithmetic
#' @export
setMethod("%%", signature(e1="ddmatrix", e2="ddmatrix"), 
  function(e1, e2){
    base.checkem(x=e1, y=e2, checks=1:3)
    if (base.ownany(dim=e1@dim, bldim=e1@bldim, ICTXT=e1@ICTXT))
      e1@Data <- e1@Data %% e2@Data
    
    return(e1)
  }
)



#' @rdname arithmetic
#' @export
setMethod("%%", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- length(e2)
    if ( (dim[1]%%len > 0 && len%%dim[1] > 0) && len > 1)
      comm.warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, ICTXT=e1@ICTXT))
      if (len==1)
        e1@Data <- e1@Data %% e2
      else {
        descx <- base.descinit(dim=e1@dim, bldim=e1@bldim, ldim=e1@ldim, ICTXT=e1@ICTXT)
        e1@Data <- base.rl2blas(x=e1@Data, descx=descx, vec=e2, FUN=5) # FUN=5 for "%%"
      }
      
    return(e1)
  }
)



#' @rdname arithmetic
#' @export
setMethod("%%", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2){
    dim <- e2@dim
    len <- length(e1)
    if ( (dim[1]%%len > 0 && len%%dim[1] > 0) && len > 1)
      comm.warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e2@bldim, ICTXT=e2@ICTXT))
      if (len==1)
        e2@Data <- e1 %% e2@Data
      else {
        descx <- base.descinit(dim=e2@dim, bldim=e2@bldim, ldim=e2@ldim, ICTXT=e2@ICTXT)
        e2@Data <- base.rl2blas(x=e2@Data, descx=descx, vec=e1, FUN=6) # FUN=0 for reverse "%%"
      }
    
    return(e2)
  }
)



#' @rdname arithmetic
#' @export
setMethod("%/%", signature(e1="ddmatrix", e2="ddmatrix"), 
  function(e1, e2){
    base.checkem(x=e1, y=e2, checks=1:3)
    if (base.ownany(dim=e1@dim, bldim=e1@bldim, ICTXT=e1@ICTXT))
      e1@Data <- e1@Data %/% e2@Data
    return(e1)
  }
)



#' @rdname arithmetic
#' @export
setMethod("%/%", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2){
    return(floor(e1 / e2))
  }
)



#' @rdname arithmetic
#' @export
setMethod("%/%", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    return(floor(e1 / e2))
  }
)
