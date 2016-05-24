### For Dirichlet distribution.

### Density of Dirichlet.
ddirichlet <- function(x, alpha, log = FALSE){
  ### length(x) should be equal to length(alpha)
  ret <- sum(lgamma(alpha)) - lgamma(sum(alpha)) + (alpha - 1) * log(x)
  if(!log){
    ret <- exp(ret)
  }
  ret
} # End of ddirichlet().

rdirichlet <- function(n, alpha){
  x <- matrix(rgamma(n * length(alpha), alpha, 1), nrow = n, byrow = TRUE)
  x / rowSums(x)
} # End of rdirichlet().
