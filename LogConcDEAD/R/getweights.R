'getweights' <- function(x) {
  if(is.vector(x)) x <- matrix(x,ncol=1)
  if(!is.matrix(x)) stop("x must be a real-valued matrix")
  xunique <- unique(x)
  w <- NULL
  for (i in 1:nrow(unique(x))) {
    w <- c(w,sum(apply(x,1,identical,xunique[i,])))
  }  
  w <- w/sum(w)
  return(list(x=xunique,w=w))
}
