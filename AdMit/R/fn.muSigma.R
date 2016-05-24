## Function which computes the location vector and the scale matrix
## __input__
## w      : [Nx1 vector] of weights
## theta  : [Nxk matrix] of draws
## __ouput__
## [list] with the following components:
## $mu    : [kx1 vector] of location
## $Sigma : [k^2x1 matrix] covariance matrix (in vector form)
## __20080429__
'fn.muSigma' <- function(w, theta, mu=NULL)
  {
    theta <- as.matrix(theta) ## univariate
    wscale <- exp(log(w)-log(sum(w)))
    if (is.null(mu))
      {
        mu <- apply(wscale*theta, 2, sum)
      }
    tmp <- theta-matrix(mu, nrow(theta), ncol(theta), byrow=TRUE)
    Sigma <- t(tmp) %*% (wscale*tmp)
    
    list(mu=as.vector(mu),
         Sigma=as.vector(Sigma))
  }
