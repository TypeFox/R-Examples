##
##  PURPOSE:   Calculate initial component allocations for a mixture model
##             or check supplied inits for consistency.
##
##             THIS IS A HELP FUNCTION, NOT TO BE CALLED BY ORDINARY USERS
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:    14/02/2010 (by taking sub-code originally included in NMixMCMC function)
##
##  FUNCTIONS:  NMixMCMCinitr
##
## ================================================================================================

## *************************************************************
## NMixMCMCinitr
## *************************************************************
##
NMixMCMCinitr <- function(z, K, w, mu, Sigma, p, n, initr)
{
    
  ##### Calculate initial r if not given
  if (missing(initr)){
    if (p == 1){
      MEANS <- matrix(rep(mu, n), ncol=K, byrow=TRUE)
      SDS   <- matrix(rep(sqrt(Sigma), n), ncol=K, byrow=TRUE)
      YY    <- matrix(rep(z, K), ncol=K)
      WW    <- matrix(rep(w, n), ncol=K, byrow=TRUE)
      PROB  <- WW * dnorm(YY, mean=MEANS, sd=SDS)
    }else{
      PROB <- matrix(0, nrow=n, ncol=K)
      for (j in 1:K){
        MEANS <- mu[((j-1)*p+1):(j*p)]
        SIGMA <- Sigma[((j-1)*p+1):(j*p),]
        PROB[,j] <- w[j] * dMVN(z, mean=MEANS, Sigma=SIGMA)        
      }        
    }
    sumPROB <- apply(PROB, 1, sum)
    sumPROB[sumPROB <= 0] <- 1
    PROB    <- PROB / matrix(rep(sumPROB, each=K), ncol=K, byrow=TRUE)
    initr <- apply(PROB, 1, which.max)
  }
  
  ##### Check and format initr
  initr <- as.numeric(initr)
  if (length(initr) != n) stop(paste("init$r must be of length ", n, sep=""))
  if (any(is.na(initr))) stop("NA in init$r")
  if (any(initr < 1) | any(initr > K)) stop(paste("init$r out of the range (must lie between ", 1, " and ", K, ")", sep=""))
  names(initr) <- paste("r", 1:n, sep="")      

  return(initr)  
}
