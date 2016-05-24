#  Get or set the coordinates (either Cartesian or fractional atoms) of an object.

coords <- function(...)
  UseMethod("coords")

'coords<-' <- function(x, value)
  UseMethod("coords<-", x)

coords.default <- function(x1, x2, x3, basis = "xyz", ...)
{
  if(!basis %in% c("xyz", "abc")) stop("Unrecognized basis")
  
  to.return <- data.frame(x1,x2,x3)
  attr(to.return, which = "basis") <- basis
  
  class(to.return) <- c("coords","data.frame")
  
  return(to.return)
}

coords.data.frame <- function(x, basis = NULL, ...)
{
  if(!is.data.frame(x)) stop("'x' must be a 'data.frame'")
  
  if(is.null(basis)){
    if(all(c("x","y","z") %in% names(x))){
      x <- x[,c("x","y","z")]
      basis <- "xyz"
    }
    else if(all(c("a","b","c") %in% names(x))){
      x <- x[,c("a","b","c")]
      basis <- "abc"
    }
    else stop("Can not convert this 'data.frame' into 'coords': Coordinates not found")
  }
  else if(!basis %in% c("xyz","abc")) stop("Unrecognized 'basis'")
  if(ncol(x) != 3L) stop("'x' must be a three-columns data.frame")

  to.return <- coords.default(x[,1], x[,2], x[,3], basis = basis)

  return(to.return)
}

coords.matrix <- function(x, basis = NULL, ...){
  if(!is.matrix(x)) stop("'x' must be a 'matrix'")
  
  to.return <- coords.data.frame(as.data.frame(x), basis = basis, ...)
  return(to.return)
}

coords.atoms <- function(x, ...)
{
  if(!is.atoms(x)) stop("'x' must be an object of class 'atoms'")
  
  to.return <- coords(x$x1,x$x2,x$x3,basis.default(x))
  return(to.return)
}

'coords<-.atoms' <- function(x, value)
{
  if(!is.atoms(x)) stop("'x' must be an object of class 'atoms'")
  
  if(!is.coords(value)) stop("'value' must be an object of class 'coords'")
  if(nrow(x) != nrow(value)) stop(paste("arguments imply different number of rows: ",nrow(x),", ",nrow(value),sep=""))
  x[c("x1","x2","x3")] <- value
  basis(x) <- basis.default(value)
  
  return(x)
}

coords.pdb <- function(x, ...)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'pdb'")
  to.return <- coords.atoms(x$atoms)
  return(to.return)
}

'coords<-.pdb' <- function(x, value)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'pdb'")
  
  if(!is.coords(value)) stop("'value' must be an object of class 'coords'")
  if(nrow(x$atoms) != nrow(value)) stop(paste("arguments imply different number of rows: ",nrow(x$atoms),", ",nrow(value),sep=""))
  coords(x$atoms) <- value
  
  return(x)
}

is.coords <- function(x)
{
  to.return <- any(class(x) == "coords")
  return(to.return)
}
