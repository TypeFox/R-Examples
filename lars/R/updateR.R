updateR <-
function(xnew, R = NULL, xold, eps = .Machine$double.eps, Gram = FALSE)
{
###Gram argument determines the nature of xnew and xold
  xtx <- if(Gram) xnew else sum(xnew^2)
  norm.xnew <- sqrt(xtx)
  if(is.null(R)) {
    R <- matrix(norm.xnew, 1, 1)
    attr(R, "rank") <- 1
    return(R)
  }
  Xtx <- if(Gram) xold else drop(t(xnew) %*% xold)
  r <- backsolvet(R, Xtx)
  rpp <- norm.xnew^2 - sum(r^2)
  rank <- attr(R, "rank")	### check if R is machine singular
  if(rpp <= eps)
    rpp <- eps
  else {
    rpp <- sqrt(rpp)
    rank <- rank + 1
  }
  R <- cbind(rbind(R, 0), c(r, rpp))
  attr(R, "rank") <- rank
  R
}

