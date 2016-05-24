llP <- function(beta, X, event, offset){
  f <- as.numeric(X%*%beta)
  ef <- exp(f)
  return(-sum(offset * ef) + sum(event * log(offset * ef)))  
}