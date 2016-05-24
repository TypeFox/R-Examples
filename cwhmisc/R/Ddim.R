Ddim <- function(x) {
  if (is.vector(x)) return( length(x) ) else  return( dim(x) )
}
