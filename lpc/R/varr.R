varr <- function(x, meanx=NULL){
  n <- ncol(x)
  p <- nrow(x)
  Y <-matrix(1,nrow=n,ncol=1)
  if(is.null(meanx)){   meanx <- rowMeans(x)}
  ans<- rep(1, p)
  xdif <- x - meanx %*% t(Y)
  ans <- (xdif^2) %*% rep(1/(n - 1), n)
  ans <- drop(ans)
  return(ans)
}
