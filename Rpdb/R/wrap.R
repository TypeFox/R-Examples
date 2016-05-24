#  Wrap atomic coordinates using periodical boundary conditions

wrap <- function(x, ...)
  UseMethod("wrap")

wrap.coords <- function(x, cryst1 = NULL, factor = NULL, ...)
{
  if(is.null(cryst1)) stop("Please specify a 'cryst1' object")
  if(is.null(factor)) factor <- 1:natom(x)
  
  if(!is.coords(x)) stop("'x' must be an object of class 'coords'")
  if(!is.cryst1(cryst1)) stop("'cryst1' must be an object of class 'cryst1'")

  b <- basis(x)
  if(b == "xyz") x <- xyz2abc(x, cryst1)
  
  centers <- centres.coords(x, factor = factor, unsplit = TRUE)
  x[centers > 1] <- x[centers > 1] - 1
  x[centers < 0] <- x[centers < 0] + 1
  
  if(b == "xyz") x <- abc2xyz.coords(x, cryst1)
  
  return(x)
}

wrap.atoms <- function(x, cryst1= NULL, factor = NULL, ...)
{
  if(is.null(cryst1)) stop("Please specify a 'cryst1' object")
  if(is.null(factor)) factor <- x$resid
  
  if(!is.atoms(x)) stop("'x' must be an object of class 'atoms'")
  if(!is.cryst1(cryst1)) stop("'cryst1' must be an object of class 'cryst1'")
  
  coords(x) <- wrap.coords(coords(x), cryst1, factor)
  
  return(x)
}

wrap.pdb <- function(x, cryst1 = x$cryst1, factor = NULL, ...)
{
  if(is.null(factor)) factor <- x$atoms$resid
  
  if(!is.pdb(x)) stop("'x' must be an object of class 'pdb'")
  if(!is.cryst1(cryst1)) stop("'cryst1' must be an object of class 'cryst1'")
  
  coords(x) <- wrap.coords(coords(x), cryst1, factor)
  
  return(x)
}
