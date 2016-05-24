#  Functions to get or set the basis of an object containing atomic coordinates.

basis <- function(x)
  UseMethod("basis")

basis.default <- function(x)
  attr(x, which = "basis")

'basis<-' <- function(x, value)
  UseMethod("basis<-", x)

'basis<-.default' <- function(x, value)
{
  if(!value %in% c("xyz","abc")) stop("Unrecognized basis")
  attr(x, "basis") <- value
  return(x)
}

basis.pdb <- function(x)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'pdb'")
  
  to.return <- basis.default(x$atoms)
  return(to.return)
}

'basis<-.pdb' <- function(x, value)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'pdb'")
  
  if(!value %in% c("xyz","abc")) stop("Unrecognized basis set")
  attr(x$atoms, "basis") <- value
  return(x)
}
