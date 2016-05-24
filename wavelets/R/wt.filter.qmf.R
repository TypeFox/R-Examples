wt.filter.qmf <- function(x, inverse=FALSE){
  L <- length(x)
  y <- x[L:1]*(-1)^(1:L) # x[L:1] is twice as fast as rev(x)
  if(inverse) y <- x[L:1]*(-1)^((1:L)-1)
  return(y)
}
