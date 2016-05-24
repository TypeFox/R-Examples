## Function which generates draws from a multivariate Student-t density
## __input__
## N     : [integer>0] number of draws
## mu    : [kx1 vector] location vector
## Sigma : [k^2x1 matrix] scale matrix (in vector form)
## df    : [integer>0] degrees of freedom
## __output__
## [Nxk matrix] of draws
## __20080427__
'fn.rmvt' <- function(N, mu, Sigma, df)
  {
    k <- length(mu)
    r <- matrix(mu, N, k, byrow=TRUE) + rmvt(N, matrix(Sigma,k,k), df)
    as.matrix(r)
  }
