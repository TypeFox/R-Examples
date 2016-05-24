##
##  PURPOSE:   EM algorithm to compute ML estimates in a normal mixture
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:                  31/12/2009
##
##  FUNCTIONS:  NMixEM
##
## ======================================================================

## *************************************************************
## NMixEM
## *************************************************************
NMixEM <- function(y, K, weight, mean, Sigma, toler=1e-5, maxiter=500)
{
  if (K < 1) stop("K must be at least 1")
  
  if (is.null(dim(y))) p <- 1
  else                 p <- ncol(y)

  if (p == 1) y <- matrix(y, ncol=1)
  else        y <- as.matrix(y)
  n <- nrow(y)

  ### Initial weights
  if (missing(weight) | K == 1) weight <- rep(1/K, K)
  if (length(weight) != K) stop(paste("weight must be of length ", K, sep=""))
  if (any(weight < 0)) stop("all weights must be non-negative")
  weight <- weight / sum(weight)

  ### Initial means
  ybar <- apply(y, 2, sum) / n
  sdy <- apply(y, 2, sd)
  if (missing(mean) | K == 1){
    step <- (6 / (K + 1)) * sdy
    mean <- matrix(ybar - 3 * sdy + step, nrow=1, ncol=p)
    if (K > 1) for (k in 2:K) mean <- rbind(mean, ybar - 3 * sdy + k * step)    
  }
  if (p == 1) mean <- matrix(mean, ncol=1)
  if (nrow(mean) != K) stop(paste("mean must have ", K, " rows", sep=""))
  if (ncol(mean) != p) stop(paste("mean must have ", p, " columns", sep=""))  

  ### Initial covariance matrix
  if (missing(Sigma) | K == 1) Sigma <- ((n-1)/n) * var(y)
  if (p == 1) Sigma <- matrix(Sigma, nrow=1, ncol=1)
  if (nrow(Sigma) != p) stop(paste("Sigma must have ", p, " rows", sep=""))
  if (ncol(Sigma) != p) stop(paste("Sigma must have ", p, " columns", sep=""))  

  ### Initial E-step, initial log-likelihood, initial value of the objective function
  fyk <- matrix(NA, nrow=n, ncol=K)
  
  for (k in 1:K) fyk[, k] <- dMVN(y, mean=mean[k,], Sigma=Sigma)
  W <- matrix(rep(weight, n), ncol=K, nrow=n, byrow=TRUE)
  wfyk <- fyk * W
  fy <- apply(wfyk, 1, sum)
  pik <- wfyk / matrix(rep(fy, K), nrow=n, ncol=K)
  pik[pik < 1e-15] <- 0
  pik <- pik / apply(pik, 1, sum)
  
  loglik <- sum(log(fy))
  
  lfyk <- log(fyk)
  lfyk[pik < 1e-15] <- -500
  Qfun   <- sum(pik * (lfyk + log(W)))
   
  ### Iterations
  iter <- 0
  not.converg <- (K > 1)
  while (iter < maxiter & not.converg){
    iter <- iter + 1
    
    ### M-step
    weight <- apply(pik, 2, mean)
    Sigma <- matrix(0, nrow=p, ncol=p)    
    for (k in 1:K){
      mean[k, ] <- apply(y * matrix(rep(pik[, k], p), ncol=p, nrow=n), 2, sum) / sum(pik[, k])
      wResid <- (y - matrix(rep(mean[k, ], n), nrow=n, ncol=p, byrow=TRUE)) * matrix(rep(sqrt(pik[, k]), p), ncol=p, nrow=n)
      Sigma <- Sigma + (t(wResid) %*% wResid) / n
    }  
    
    ### E-step, evaluate log-likelihood and the objective function
    for (k in 1:K) fyk[, k] <- dMVN(y, mean=mean[k,], Sigma=Sigma)
    W <- matrix(rep(weight, n), ncol=K, nrow=n, byrow=TRUE)
    wfyk <- fyk * W
    fy <- apply(wfyk, 1, sum)
    pik <- wfyk / matrix(rep(fy, K), nrow=n, ncol=K)
    pik[pik < 1e-15] <- 0
    pik <- pik / apply(pik, 1, sum)
    
    loglik <- c(loglik, sum(log(fy)))

    lfyk <- log(fyk)
    lfyk[pik < 1e-15] <- -500    
    Qfun   <- c(Qfun, sum(pik * (lfyk + log(W))))
    
    #criterion <- abs((loglik[iter + 1] - loglik[iter]) / loglik[iter + 1])
    criterion <- abs((Qfun[iter + 1] - Qfun[iter]) / Qfun[iter + 1])     

    if (criterion < toler) not.converg <- FALSE
  }  

  if (p == 1){
    mean <- as.numeric(mean)
    Sigma <- as.numeric(Sigma)
  }  

  ### AIC, BIC
  nPar <- (K - 1) + K * p + p * (1 + p) / 2   ## number of free parameters
  aic <- -2 * loglik[length(loglik)] + 2 * nPar
  bic <- -2 * loglik[length(loglik)] + log(n) * nPar
  
  RET <- list(K           = K,
              weight      = weight,
              mean        = mean,
              Sigma       = Sigma,
              loglik      = loglik[length(loglik)],
              aic         = aic,
              bic         = bic,
              iter        = iter,
              iter.loglik = loglik,
              iter.Qfun   = Qfun,
              dim         = p,
              nobs        = n)

  class(RET) <- "NMixEM"
  return(RET)  
}  
