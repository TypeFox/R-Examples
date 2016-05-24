as.array3 <- function(x){
#  as.array3 produces a 3-d array from a vector, matrix or array of 
#  up to 3 dimensions, preserving names.

  dimx <- dim(x)
  ndim <- length(dimx)
#   If dimension is already 3
  if (ndim==3) return(x)
#   Otherwise, set up an error message dimension higher than 3
  xName <- substring(deparse(substitute(x)), 1, 22) 
  if (ndim>3)
    stop('length(dim(', xName, ") = ", ndim, ' > 3')
#  If dimension less than 3, ...
  x.      <- as.matrix(x)  #  coerce to matrix
  xNames  <- dimnames(x.)  #  get dimension names if any
  dim(x.) <- c(dim(x.), 1) #  add a unit value for 3rd dimension
#  Assign dimension names
  if(is.list(xNames))
    dimnames(x.) <- list(xNames[[1]], xNames[[2]], NULL)
#  Return result
  x. 
}

