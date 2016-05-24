gcv.genlasso <- function(object) {
  lams = object$lambda
  df = object$df
  n = length(object$y)
  ymat = matrix(object$y,n,length(lams))
  pred = object$fit

  err = colMeans((ymat-pred)^2)/(1-df/n)^2
  names(err) = round(lams,3)
  lam.min = lams[which.min(err)]
  
  out = list(err=err,lambda=lams,lambda.min=lam.min)
  class(out) = c("gcv.genlasso", "list")
    
  return(out)  
}
