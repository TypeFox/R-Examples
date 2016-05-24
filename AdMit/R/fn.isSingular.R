## Function which checks whether a square matrix is singular
## __input__
## A   : [kxk matrix]
## tol : [double] tolerance for the determinant (default: 1e15)
## __output__
## [logical] is the matrix singular?
## __20080502__
'fn.isSingular' <- function(A, tol=1e25)
  {
    tmp <- abs(det(A))
    as.logical(tmp>=tol | tmp<=1/tol)
  }
    
