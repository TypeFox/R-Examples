Predict.matrix.nu.smooth <- function(object, data)
{ 
  X <- Predict.matrix.pspline.smooth(object,data)
  if(is.null(object$xt$weights))
    h <- rep(1, nrow(X))
  else
    h <- object$xt$weights
  R <- NULL
  for(i in 1L:ncol(X)) {
    x <- X[,i] * h
    R <- cbind(R, x / sum(x))
  }
	
  return(R)
}

