"rwishart" <-
function(N, df, Sigma)
  { p = nrow(Sigma)
    SqrtSigma <- t(chol(Sigma))
    tmp <- array(0, c(p,p,N))
    for (i in 1:N)
      {
        Z <- matrix(rnorm(df*p),nrow=p,ncol=df)
        ZS <- crossprod(Z,SqrtSigma)
        tmp[,,i] <- crossprod(ZS)
      }
    if(N==1)
      {
        return(matrix(tmp,p,p))
      }
    else
      {
        return(tmp)
      }
  }

# log density for a Wishart
# Based on code from MCMCpack

"ldwishart" <- function(W, v, S)
{
    k <- nrow(S)
    lgammapart <- 0
    for (i in 1:k) {
        lgammapart <- lgammapart + lgamma((v + 1 - i)/2)
    }
    denom <- lgammapart + (v * k/2)*log(2) + (k * (k - 1)/4)*log(pi)
    detS <- determinant(S)$modulus[1]
    detW <- determinant(W)$modulus[1]
    hold <- solve(S) %*% W
    tracehold <- sum(hold[row(hold) == col(hold)])
    num <- detS*(-v/2) + detW*((v - k - 1)/2) + (-1/2 * tracehold)
    return(num-denom)
}
