rkwp <- function(n, k, beta, a, b)
{
  u <- runif(n, 0, 1)
  return(beta/((1-(1-(1-u)^(1/b))^(1/a))^(1/k)))
}