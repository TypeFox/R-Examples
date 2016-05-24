#  Create an object of class 'pdb'

pdb <- function(...)
  UseMethod("pdb")

pdb.default <- function(atoms, cryst1 = NULL, conect = NULL, remark = NULL, title = NULL, ...)
{
  if(missing(atoms)) stop("Please specify at least a 'atoms' object")
  if(!is.atoms(atoms)) stop("'atoms' must be an object of class 'atoms'")
  
  if(!is.null(cryst1) & !is.cryst1(cryst1)) stop("'cryst1' must be an object of class 'cryst1'")
  if(!is.null(conect) & !is.conect(conect)) stop("'conect' must be an object of class 'conect'")
  
  if(is.list(title ) | !is.null(dim(title ))) stop("'title' must be a vector of character strings")
  if(is.list(remark) | !is.null(dim(remark))) stop("'remark' must be a vector of character strings")
  
  if(!is.character(title ) & !is.null(title )) title  <- as.character(title )
  if(!is.character(remark) & !is.null(remark)) remark <- as.character(remark)
  
  to.return <- list(title = title, remark = remark, cryst1 = cryst1, atoms = atoms, conect = conect)
  
  class(to.return) <- c("pdb","list")
  return(to.return)
}

is.pdb <- function(x)
{
  to.return <- any(class(x) == "pdb")
  return(to.return)
}
