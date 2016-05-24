outlier.prob.metaplus <- function(object) {
  isreg <- !is.null(object$mods)

  muhat <- object$results[1,1]
  tau2 <- object$results[2,1] 
  tau2out <- object$results[3,1] 
  lpoutlier <- object$results[4,1]
  if (isreg) xcoef <- object$results[5:dim(object$results)[1],1]
                          
  w <- 1.0/(tau2+object$sei^2)
  if (isreg) ll1 <- -0.5*(log(2*pi)-log(w)+w*(object$yi-muhat-as.vector(object$mods %*% matrix(xcoef,ncol=1)))^2)
  else ll1 <- -0.5*(log(2*pi)-log(w)+w*(object$yi-muhat)^2)
  w <- 1.0/(tau2out+object$sei^2)
  if (isreg) ll2 <- -0.5*(log(2*pi)-log(w)+w*(object$yi-muhat-as.vector(object$mods %*% matrix(xcoef,ncol=1)))^2)
  else ll2 <- -0.5*(log(2*pi)-log(w)+w*(object$yi-muhat)^2)
  poutlier <- exp(lpoutlier)/(1+exp(lpoutlier))
  p <- c(1-poutlier,poutlier)
  ll <- cbind(ll1,ll2)
  l <- exp(ll)

  prop <- t(p*t(l))/apply(t(p*t(l)),1,sum)
  return(list(outlier.prob=prop[,2]))
}