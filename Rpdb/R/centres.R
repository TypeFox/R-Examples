#  Compute geometric and mass centres of groups of atoms.

centres <- function(...)
  UseMethod("centres")

centres.coords <- function(x, factor = NULL, weights = NULL, unsplit = FALSE, na.rm = FALSE, ...)
{
  if(!is.coords(x)) stop("'x' must be an object of class 'coords'")

  if(is.null(factor )) factor  <- rep("",natom(x))
  if(is.null(weights)) weights <- rep(1 ,natom(x))

  w  <- split(weights, factor)
  x1 <- split(x$x1   , factor)
  x2 <- split(x$x2   , factor)
  x3 <- split(x$x3   , factor)
  
  w.mean <- function(x, w)
    sum(x*w/sum(w, na.rm = na.rm))
  
  x1.mean <- mapply(w.mean, x1, w, SIMPLIFY = FALSE)
  x2.mean <- mapply(w.mean, x2, w, SIMPLIFY = FALSE)
  x3.mean <- mapply(w.mean, x3, w, SIMPLIFY = FALSE)
  
  if(unsplit){
    x1.mean <- unsplit(x1.mean, factor)
    x2.mean <- unsplit(x2.mean, factor)
    x3.mean <- unsplit(x3.mean, factor)
  }
  
  x1.mean <- unlist(x1.mean)
  x2.mean <- unlist(x2.mean)
  x3.mean <- unlist(x3.mean)
  
  to.return <- coords.default(x1.mean, x2.mean, x3.mean, basis(x))
  
  return(to.return)
}

centres.atoms <- function(x, factor = NULL, weights = NULL, unsplit = FALSE, na.rm = FALSE, ...)
{
  if(!is.atoms(x)) stop("'x' must be an object of class 'atoms'")
  
  if(is.null(factor)) factor <- x$resid
  
  to.return <- centres.coords(coords(x), factor, weights, unsplit, na.rm)
  return(to.return)
}

centres.pdb <- function(x, factor = NULL, weights = NULL, unsplit = FALSE, na.rm = FALSE, ...)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'pdb'")
  
  to.return <- centres.atoms(x$atoms, factor, weights, unsplit, na.rm)
  return(to.return)
}
