smooth.construct.nu.smooth.spec <-
function(object, data, knots)
{
  if(is.null(object$xt$weights))
    h <- NULL
  else
    h <- object$xt$weights
  B <- smooth.construct.ps.smooth.spec(object, data, knots)
  class(B) <- "nu.smooth"
  X <- B$X
  if(is.null(h))
    h <- rep(1, nrow(X))
  R <- NULL
  for(i in 1L:ncol(X))  {
    x <- X[,i] * h
    R <- cbind(R, x/sum(x))
  }
  B$X <- R

  return(B)
}

