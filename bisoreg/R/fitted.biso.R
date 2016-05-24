`fitted.biso` <- function(object,xnew,...){
  ## PURPOSE:
  ## Compute fitted values for a "biso" object
  ##
  ## INPUTS:
  ## obj  = a "biso" object
  ## xnew = a vector of new x-values at which to
  ##        to compute predicted values
  ##
  ## OUTPUTS:
  ## yhat = vector of predicted values
  ##
  o <- order(object$x)
  ro <- sapply(1:length(object$x), function(x) which(o==x))
  x <- object$x[o]
  y <- object$y[o]
  m <- object$m
  u <- as.matrix(object$postdraws[,1:m])
  u0 <- object$postdraws$u0
  if(missing(xnew)) xnew <- x
  wnew <- get.W((xnew-object$xmin)/object$xrng,m)
  yhat <- sapply(1:nrow(wnew),function(i) median(u0 + u%*%wnew[i,]))
  return(yhat[ro])
  }

