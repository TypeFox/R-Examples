"as.xyz" <- function(x) {
  if(is.vector(x))
    x = matrix(x, nrow=1)
  
  dims <- dim(x)
  if(!(dims[2L]%%3==0))
    warning("number of cartesian coordinates not a multiple of 3")
  
  class(x) <- c("xyz", "matrix")
  return(x)
}

