translate <- function(obj)
{
  gamma <- coef(obj)
  sigma <- obj$scale
  beta <- -gamma/sigma
  scale <- sigma
  list(beta = beta, scale = 1/sigma)
}