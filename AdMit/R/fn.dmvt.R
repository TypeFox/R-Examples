## Function which computes the density values of a multivariate Student-t density
## __input__
## x     : [Nxk matrix] of values
## mu    : [kx1 vector] of mean
## Sigma : [k^2x1 matrix] scale matrix (in vector form)
## df    : [integer>0] degrees of freedom parameter
## log    : [logical] natural logarithm output 
## __output__
## [Nx1 vector] of density values
## __20080429__
'fn.dmvt' <- function(x, mu, Sigma, df, log)
  {
    k <- length(mu)
    r <- dmvt(as.matrix(x), as.vector(mu), matrix(Sigma,k,k), df, log)
    as.vector(r)
  }
