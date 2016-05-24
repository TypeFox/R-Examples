rmnorm <-
function(n, mean=NULL, vcov=1)
{
  d <- NCOL(vcov)
  y <- matrix(rnorm(n * d), d, n)
  y <- crossprod(y,chol(vcov))
  if(!is.null(mean)){ y <- t(mean + t(y)) }
  return(y)
}
