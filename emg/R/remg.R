remg <- function(n, mu=0, sigma=1, lambda=1)
{
  rnorm(n, mu, sigma) + rexp(n, lambda)
}