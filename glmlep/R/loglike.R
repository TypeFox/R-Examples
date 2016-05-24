loglike <-
function(x,y,beta,family=c("gaussian","binomial")){
  
  ## -2 log-lieklihood for different family
  
  family = match.arg(family)
  if(any(x[,1]!=1)) x <- cbind(1,x)
  n <- length(y)
  loss <- switch(family, gaussian = n*log(sum((y - x%*%beta)^2)), binomial = -2*(crossprod(y,x%*%beta) - sum(log(1+exp(x%*%beta)))))
  return(loss)  
}
