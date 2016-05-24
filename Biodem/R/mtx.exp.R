"mtx.exp" <- function(X,n){
## Function to calculate the n-th power of a matrix X
  if(n != round(n)) {
    n <- round(n)
    warning("rounding exponent `n' to", n)
  }
  phi <- diag(nrow = nrow(X))
  pot <- X # the first power of the matrix.
  while (n > 0)
    {
      if (n %% 2)
        phi <- phi %*% pot
      n <- n %/% 2
      pot <- pot %*% pot
    }
  return(phi)
}


 
