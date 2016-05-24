Rsquared <- function(object, pi.hat)
{
# objects of class ppar
# computes Pearson R^2 and SS R^2 for objects of class ppar

  #Pi <- pmat(object)                              #expected values
  if (length(object$pers.ex) > 0){
    y <- as.vector(t(object$X[-object$pers.ex,])) #observed values
  } else {
    y <- as.vector(t(object$X))
  }
  pi.hat <- as.vector(t(pi.hat))

  R.P <- cor(y, pi.hat)^2                                 #Squared Pearson correlation
  R.SS <- 1-(sum((y - pi.hat)^2)/sum((y - mean(y))^2))    #SS-R^2
  
  loglik.full <- sum(y*log(pi.hat) + (1-y)*log(1-pi.hat), na.rm = TRUE)  #full likelihood
  loglik.0 <- sum(y*log(mean(y)) + (1-y)*log(1-mean(y)))    #L0 (Agresti, Sec. 6.2.5)
  R.MF <- (loglik.0 - loglik.full)/loglik.full

  result <- list(R2.P = R.P, R2.SS = R.SS, R2.MF = R.MF)
  result
}