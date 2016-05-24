predict.composition <- function(object, newdata, n.method = "median", l.bound = NULL, ...){
  
  x <- newdata
  X <- model.matrix(object$n.formula, x)
  b <- object$b
  v <- matrix(0, nrow(b), object$D + 1)
  K <- ncol(X)
  n <- nrow(X)
  if (is.null(l.bound)) l.bound <- object$l.bound
  l.bound <- ifelse(length(l.bound)==1, rep(l.bound, n), l.bound)
  
  if(l.bound > object$D + 1) stop("Lower bound not valid.")
  
  for (i in 1:nrow(b)){
    
    if (n.method=="median") {
      n <- min(which(cumsum(dtnegbin(0:(object$D+10), mu=exp(X%*%b[i, 1:K]), dispersion=exp(b[i, K+1]), l.bound = l.bound)) >= .5)[1], object$D + 1)
    }
    
    if (n.method=="mode") {
      n <- which.max(dtnegbin(0:(object$D + 1), mu=exp(X%*%b[i, 1:K]), dispersion=exp(b[i, K+1]), l.bound = l.bound))
    }
    
    G <- model.matrix(object$v.formula, data.frame(x, n=n))
    y  <- G%*%object$g[,,i]
    v[i, 1:n]  <- y.v(y[1:(n-1)])
    
  }
  
  class(v) <- "comphat"
  v
}
