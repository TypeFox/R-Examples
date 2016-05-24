## partial derivatives of log complete loglikelihood but using the multivariate
## logit transform for the PIs.

### partial.q. Return a matrix with dimension M * N.
partial.q <- function(x, PI, MU, S, t, logit.PI = TRUE){
  ##
  ## t is the posterior classification probability.
  ## PI is the mixing component proportions
  ## MU, S are means and variance-covariance matrices of the mixture components
  ## 

  K <- length(PI)
  N <- dim(x)[1]
  p <- dim(x)[2]

  mp <- p * (p + 1) / 2
  M <- K - 1 + K * p + K * p * (p + 1) / 2
  G <- Gmat(p)
  Id <- diag(p)

  invS <- apply(S, 3, solve)
  dim(invS) <- c(p, p, K)

  SS <- lapply(1:N, function(i){
    Sbuf.1 <- NULL
    if(K != 1){
      if(logit.PI){    # partial wrt the multivariate logit parameter
        Sbuf.1 <- t[i, -K] - PI[-K]
      } else{
        tmp <- t[i, ] / PI
        Sbuf.1 <- tmp[1:(K-1)] - tmp[K]
      }
    }

    Sbuf.2 <- list()
    Sbuf.3 <- list()
    for(j in 1:K){
      tmp <- x[i,] - MU[j,]
      dim(tmp) <- c(p, 1)
      tmp.1 <- t[i,j] * invS[,,j]
      Sbuf.2[[j]] <- tmp.1 %*% tmp
      Smat <-  tmp.1 %*% (tmp %*% t(tmp) %*% invS[,,j] - Id) / 2
      Sbuf.3[[j]] <- as.vector(Smat) %*% G
    }

    c(Sbuf.1, do.call("c", Sbuf.2), do.call("c", Sbuf.3))
  })

  do.call("cbind", SS)
} # End of partial.q().

