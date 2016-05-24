tritomat <-
function(dmat,n)
{ ## dmat : lower trangle elements
  ## n    : dimension of the square matrix
  ## value : matrix with diagonal zero filled and 
  ##       upper triangle mapped
  out <- matrix(0,n,n)
  out[lower.tri(out)] <- dmat
  out <- t(out)
  out[lower.tri(out)] <- dmat
  t(out)
}
