#
# Simulate multivariate normal variable. Each ROW is a random vector

.rmvn <- function(n, mu, sigma, isChol = FALSE)
{
  d <- length(mu)
  
  if(isChol) cholFact <- t(sigma) else cholFact <- t( chol(sigma) )

  return( t( mu + cholFact %*% matrix(rnorm(n*d), d, n) ) ) 
  
}
  