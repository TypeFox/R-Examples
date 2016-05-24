SetLambda <-
function(x, y, lambda.min, nlambda, penalty.factor){
  
  ## generate lambda sequence
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  ind <- which(penalty.factor!=0)
  if (length(ind)!=p) {
    fit <- lm(y~x[, -ind])
  } else fit <- lm(y~1)
  
  z <- crossprod(x[,ind, drop=FALSE], fit$residuals) / n
  
  lambda.max <- max(abs(z)/(penalty.factor[ind]))
  
  if (lambda.min==0) {
    lambda <- c(exp(seq(log(lambda.max),log(.001*lambda.max),len=nlambda-1)),0)
  } else lambda <- exp(seq(log(lambda.max),log(lambda.min*lambda.max),len=nlambda))
  
  return(lambda)
}
