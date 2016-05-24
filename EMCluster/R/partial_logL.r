## functions to calculate partial derivatives of loglikelihood function of the
## multivariate Gaussian mixture model but using the multivariate logit
## transform for the PIs.

### partial.logL. Return a matrix with dimension M * N.
partial.logL <- function(x, PI, MU, S, mixdensvals = NULL, logit.PI = TRUE){
  ##
  ## PI is the vector of mixing component proportions
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
  LTS <- variance2LTSigma(S)

  if(is.null(mixdensvals)){
    mixdensvals <- dmixmvn(x = x, pi = PI, Mu = MU, LTSigma = LTS)
  }

  SS <- lapply(1:N, function(i){
    mixdens <- mixdensvals[i]
    gauscomps <- unlist(lapply(1:K, function(k){
                                      dmvn(x = x[i, ], mu = MU[k, ],
                                           LTsigma = LTS[k, ]) }))

    Sbuf.1 <- NULL
    if(K != 1){
      if(logit.PI){    # partial wrt the multivariate logit parameter
        Sbuf.1 <- PI[-K] * gauscomps[-K] / mixdens - PI[-K]
      } else{
        Sbuf.1 <- (gauscomps[-K] - gauscomps[K]) / mixdens
      }
    }

    Sbuf.2 <- list()
    Sbuf.3 <- list()
    for(j in 1:K){
      tmp <- x[i,] - MU[j,]
      dim(tmp) <- c(p, 1)
      Sbuf.2[[j]] <- -1 * invS[,,j] %*% tmp * PI[j] * gauscomps[j] / mixdens
      Smat <- 0.5 * invS[,,j] %*% (tmp %*% t(tmp) %*% invS[,,j] - Id) *
              gauscomps[j] / mixdens
      Sbuf.3[[j]] <- as.vector(Smat) %*% G
    }

    c(Sbuf.1, do.call("c", Sbuf.2), do.call("c", Sbuf.3))
  })

  do.call("cbind", SS)
} # End of partial.logL().

