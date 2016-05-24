#  Evaluate the number of atom in an object containing atomic coordinates.

natom <- function(x, ...)
  UseMethod("natom")

natom.coords <- function(x, factor = NULL, ...)
{
  if(length(factor) == 0)
    nrow(x)
  else
    unlist(lapply(split(x, factor), nrow))
}

natom.atoms <- function(x, factor = NULL, ATOM = TRUE, HETATM = TRUE, ...)
{
  M <- ATOM & x$recname == "ATOM" | HETATM & x$recname == "HETATM"
  x <- x[M,]
  if(length(factor) != 0) factor <- factor[M]
  NextMethod("natom", x)
}

natom.pdb <- function(x, factor = NULL, ATOM = TRUE, HETATM = TRUE, ...)
  natom.atoms(x$atoms, factor, ATOM, HETATM, ...)