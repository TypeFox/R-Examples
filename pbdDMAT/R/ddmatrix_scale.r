# this is a complete travesty, but I'm not sure how to make it simpler, really
# hope you like spaghetti

### Utility functions for scale
dmat.clmn <- function(x, na.rm=TRUE)
{
  if (x@dim[1L]==1)
    return(x@Data)
  
  Data <- colSums(x@Data, na.rm=na.rm) / as.double(x@dim[1L])
  dim(Data) <- c(1L, base::length(Data))
  
  # local dimension correction for reduction
  ldim <- dim(Data)
  ldim[2L] <- dmat.allcolreduce(ldim[2L], op='max', ICTXT=x@ICTXT)
  
  if (dim(Data)[2L] != ldim[2L])
    Data <- matrix(0.0, ldim[1L], ldim[2L])
  
  out <- dmat.allcolreduce(Data, op='sum', ICTXT=x@ICTXT)
  
  return( out )
}


dmat.clscl <- function(x, na.rm=TRUE)
{
  if (x@dim[1L]==1)
    return(abs(x@Data))
  
  n <- x@ldim[2L]
  
  len <- integer(n)
  Data <- matrix(0.0, nrow=1L, ncol=n)
  for (i in 1L:x@ldim[2L]){
    v <- x@Data[, i]
    v <- v[!is.na(v)]
    len[i] <- length(v)
    Data[1L, i] <- sum(v^2)
  }
  
  # local dimension correction for reduction
  maxn <- dmat.allcolreduce(n, op='max', ICTXT=x@ICTXT)
  
  if (dim(Data)[2L] != maxn)
    Data <- matrix(0.0, 1L, maxn)
  
  if (!base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
    len <- integer(maxn)
  len <- dmat.allcolreduce(len, op='sum', ICTXT=x@ICTXT)
  
  len <- sapply(len, function(i) max(i-1L, 1L))
  
  Data <- base::sweep(Data, len, FUN="/", MARGIN=2, check.margin=FALSE)
  out <- sqrt( dmat.allcolreduce(Data, op='sum', ICTXT=x@ICTXT) )
  
  return( out )
}


### logical
dmat.scale.center.logical <- function(x)
{
  cntr <- dmat.clmn(x, na.rm=FALSE)
  iown <- base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT)
  
  if (iown)
    x@Data <- base::scale(x@Data, center=cntr, scale=FALSE)
  
  # attributes
  dim <- c(1, x@dim[2L])
  ldim <- base.numroc(dim=dim, bldim=x@bldim, ICTXT=x@ICTXT)
  if (iown){
    Data <- cntr
    dim(Data) <- c(1L, base::length(cntr))
  }
  else
    Data <- matrix(0.0, 1, 1)
  
  center <- new("ddmatrix", Data=Data, dim=dim, ldim=ldim, bldim=x@bldim, ICTXT=x@ICTXT)
  attr(x@Data, "scaled:center") <- center
  
  return( x )
}

dmat.scale.scale.logical <- function(x)
{
  scl <- dmat.clscl(x, na.rm=FALSE)
  
  iown <- base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT)
  
  if (iown)
    x@Data <- base::scale(x@Data, center=FALSE, scale=scl)
  else {}
  
  # attributes
  dim <- c(1, x@dim[2L])
  ldim <- base.numroc(dim=dim, bldim=x@bldim, ICTXT=x@ICTXT)
  if (iown){
    Data <- scl
    dim(Data) <- c(1L, base::length(Data))
  }
  else
    Data <- matrix(0.0, 1, 1)
  
  scale <- new("ddmatrix", Data=Data, dim=dim, ldim=ldim, bldim=x@bldim, ICTXT=x@ICTXT)
  attr(x@Data, "scaled:scale") <- scale
  
  return( x )
}

### matrix/vector
dmat.scale.center.atomic <- function(x, center)
{
  if (is.matrix(center))
    dim(center) <- NULL
  
  descx <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
  ret <- base.pdsweep(x=x@Data, descx=descx, vec=center, MARGIN=2, FUN="-")
  x@Data <- ret
  attr(x@Data, "scaled:center") <- center
  
  return( x )
}

dmat.scale.scale.atomic <- function(x, scale)
{
  if (is.matrix(vector))
    dim(scale) <- NULL
  
  descx <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
  ret <- base.pdsweep(x=x@Data, descx=descx, vec=scale, MARGIN=2L, FUN="-")
  x@Data <- ret
  attr(x@Data, "scaled:scale") <- scale
  
  return( x )
}


### ddmatrix
dmat.scale.center.ddmatrix <- function(x, center)
{
  if (prod(center@dim) != x@dim[2L])
    comm.stop("length of 'center' must equal the number of columns of 'x'")
  else if (any(center@bldim != x@bldim))
    comm.stop("distributed matrices 'x' and 'center' must have the same blocking dimension")
  else if (center@ICTXT != x@ICTXT)
    comm.stop("distributed matrices 'x' and 'center' must belong to the same ICTXT")
  else if (center@dim[1L] != 1 && center@dim[2L] != 1)
    comm.stop("can't do this yet") #FIXME
  
  if (center@dim[2L] == 1)
    center <- t(center)
  
  ldim <- center@ldim
  ldim <- dmat.allcolreduce(ldim, op='max', ICTXT=center@ICTXT)
  
  if (ownany(dim=center@dim, bldim=center@bldim, ICTXT=center@ICTXT))
    cntr <- center@Data
  else
    cntr <- matrix(0.0, ldim[1L], ldim[2L])
  
  cntr <- dmat.allcolreduce(cntr, op='sum', ICTXT=x@ICTXT)
  
  x@Data <- base::scale(x@Data, center=cntr, scale=FALSE)
  attr(x@Data, "scaled:center") <- center
  
  return( x )
}

dmat.scale.scale.ddmatrix <- function(x, scale)
{
  if (prod(scale@dim) != x@dim[2L])
    comm.stop("length of 'scale' must equal the number of columns of 'x'")
  else if (any(scale@bldim != x@bldim))
    comm.stop("distributed matrices 'x' and 'scale' must have the same blocking dimension")
  else if (scale@ICTXT != x@ICTXT)
    comm.stop("distributed matrices 'x' and 'scale' must belong to the same BLACS context")
  else if (scale@dim[1L] != 1 && scale@dim[2L] != 1)
    comm.stop("can't do this yet") #FIXME
  
  if (scale@dim[2L] == 1)
    scale <- t(scale)
  
  ldim <- scale@ldim
  ldim <- dmat.allcolreduce(ldim, op='max', ICTXT=scale@ICTXT)
  
  if (ownany(dim=scale@dim, bldim=scale@bldim, ICTXT=scale@ICTXT))
    scl <- scale@Data
  else
    scl <- matrix(0.0, ldim[1L], ldim[2L])
  
  scl <- dmat.allcolreduce(scl, op='sum', ICTXT=x@ICTXT)
  
  x@Data <- base::scale(x@Data, center=FALSE, scale=scl)
  attr(x@Data, "scaled:scale") <- scale
  
  return( x )
}



# Abandan hope all ye who enter in


#' Scale
#' 
#' Centers and/or scales the columns of a distributed matrix.
#' 
#' @param x 
#' numeric distributed matrix.
#' @param center 
#' logical value, determines whether or not columns are zero
#' centered
#' @param scale 
#' logical value, determines whether or not columns are rescaled
#' to unit variance
#' 
#' @return 
#' Returns a distributed matrix.
#' 
#' @keywords Methods
#' 
#' @examples
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' comm.set.seed(diff=T)
#' 
#' x <- ddmatrix("rnorm", 10, 10)
#' y <- scale(x)
#' 
#' print(y)
#' 
#' finalize()
#' }
#' 
#' @name ddmatrix-scale
#' @rdname ddmatrix-scale
NULL



#' @rdname ddmatrix-scale
#' @export
setGeneric(name = "scale", useAsDefault = base::scale, package="pbdDMAT")



#' @rdname ddmatrix-scale
#' @export
setMethod("scale", signature(x="ddmatrix", center="ANY", scale="ANY"),
  function(x, center=TRUE, scale=TRUE)
  {
    ### Center
    # ddmatrix
    if (is.ddmatrix(center))
      x <- dmat.scale.center.ddmatrix(x=x, center=center)
    # logical
    else if (is.logical(center)){
      if (center)
        x <- dmat.scale.center.logical(x=x)
    }
    # global matrix/vector
    else if (is.matrix(center) || (is.vector(center))){
      x <- dmat.scale.center.atomic(x=x, center=center)
    }
    # error
    else 
      comm.stop("ERROR : invalid argument for 'center'")
    
    ### Scale
    # ddmatrix
    if (is.ddmatrix(scale))
      x <- dmat.scale.scale.ddmatrix(x=x, scale=scale)
    # logical
    else if (is.logical(scale)){
      if (scale)
        x <- dmat.scale.scale.logical(x=x)
    }
    # global matrix/vector
    else if (is.matrix(scale) || (is.vector(scale) ))
      x <- dmat.scale.scale.atomic(x=x, scale=scale)
    # error
    else 
      comm.stop("ERROR : invalid argument for 'scale'")
    
    return( x )
  }
)



