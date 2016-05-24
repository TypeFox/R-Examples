`af` <- function(x) {
  if(is.vector(x)) {
   p <- (x[1]+0.5*x[2])/sum(x)
  } else
  if(is.matrix(x)) {
   p <- (x[,1]+0.5*x[,2])/apply(x,1,sum)
  }
  return(p)
}

