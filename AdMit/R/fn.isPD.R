## Function which checks whether a square matrix is positive definite
## __input__
## A   : [kxk matrix]
## __output__
## [logical] is the matrix positive definite?
## __20080427__
'fn.isPD' <- function(A)
  {
    as.logical(all(eigen(A)$values>0))
  }
    
