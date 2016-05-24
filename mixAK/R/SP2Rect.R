##
##  PURPOSE:   Create a symmetric matrix from its lower triangle
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   16/01/2008
##
##  FUNCTION:  SP2Rect
##
## ======================================================================

## *************************************************************
## SP2Rect
## *************************************************************
SP2Rect <- function(LT, dim)
{
  dim <- dim[1]
  if (dim <= 0) stop("p must be positive")
  
  lLT <- (dim*(dim + 1))/2
  if (length(LT) != lLT) stop(paste("LT must be of length ", lLT, sep=""))
    
  A <- diag(dim)
  A[lower.tri(A, diag=TRUE)] <- LT
  A[upper.tri(A, diag=FALSE)] <- t(A)[upper.tri(A, diag=FALSE)]

  return(A)    
}
