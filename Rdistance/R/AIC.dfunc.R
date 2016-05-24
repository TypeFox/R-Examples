AIC.dfunc=function (object, ..., k = 2, n = length(object$dist)) 
{
  p <- length(coef(object))
  k*p + 2*object$loglik + (k * p * (p + 1))/(n - p - 1)
}