# funcao para simular 
rNCDun <- function(N = 1, nu = Inf, rho, delta, two.sided=TRUE)
{
  R <- rho^0.5%*%t(rho^0.5)+diag(1-rho)
  maxabs <- function(x) return(max(abs(x)))
  if (two.sided==TRUE)
  {   
    if (nu==Inf) x <- apply(rmultvariate(N, delta, R), 1, maxabs) else
      x <- apply(rmultvariate(N, delta, R), 1, maxabs)/sqrt(rchisq(N, nu)/nu)
  } else
  {
    if (nu==Inf) x <- apply(rmultvariate(N, delta, R), 1, max) else
      x <- apply(rmultvariate(N, delta, R), 1, max)/sqrt(rchisq(N, nu)/nu)
  }   
  return(x)
}

