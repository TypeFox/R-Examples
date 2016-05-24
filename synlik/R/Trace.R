# Finds the trace of a _SQUARE_ matrix
.Trace <- function(m) 
{
  if ( !is.matrix(m) || (dim(m)[1] != dim(m)[2]) ) stop("m must be a square matrix")
  
  return( sum(diag(m)) )
}