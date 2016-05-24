cwdeviance <- function(object, pi.hat)
{
# computes casewise deviance for objects of class ppar

  X <- object$X.ex
  loglik.full <- sum(X*log(pi.hat)+(1-X)*log(1-pi.hat), na.rm = TRUE)  #for ordinary logistic regression
  npar.full <- (dim(object$W)[2])+sum(object$npar)           #number of estimated item + person parameters
  npar.sat <- sum(nrow(pi.hat)*ncol(pi.hat))

  value <- -2*loglik.full
  df <- npar.sat-npar.full
  p.value <- 1-pchisq(value, df = df)
  
  result <- list(value = value, df = df, p.value = p.value)
  result
}