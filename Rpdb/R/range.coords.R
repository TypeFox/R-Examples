#  Determine the range of atomic coordinates.

range.coords <- function(x, na.rm = FALSE, finite = FALSE, ...)
{
  if(!is.coords(x)) stop("'x' must be an object of class 'coords'")

  to.return <- lapply(x[,c("x1","x2","x3")], range, na.rm, finite)
  to.return <- as.data.frame(to.return, row.names = c("min","max"))

  c.names <- unlist(strsplit(basis(x), ""))
  colnames(to.return) <- c.names
  
  return(to.return)
}

range.atoms <- function(x, na.rm = FALSE, finite = FALSE, ...)
{
  if(!is.atoms(x)) stop("'x' must be an object of class 'atoms'")
  
  to.return <- range.coords(coords(x), na.rm, finite)
  return(to.return)
}

range.pdb <- function(x, na.rm = FALSE, finite = FALSE, ...)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'pdb'")
  
  to.return <- range.atoms(x$atoms, na.rm, finite)
  return(to.return)
}