# Sampling from Truncated multivariate t distribution using 
#
# a) Rejection sampling
# b) Gibbs sampling
# 
# Author: Stefan Wilhelm, Manjunath B G
#
# Literatur:
# (1) Rejection Sampling : None
# (2) Gibbs Sampling : 
# Geweke (1991) "Efficient simulation from the multivariate normal and Student-t distributions 
# subject to linear constraints and the evaluation of constraint probabilities"
###############################################################################

rtmvt <- function(n, mean = rep(0, nrow(sigma)), 
  sigma = diag(length(mean)), 
  df = 1, 
  lower = rep(-Inf, length = length(mean)), 
  upper = rep( Inf, length = length(mean)),
  algorithm=c("rejection", "gibbs"), ...) {
  
  algorithm <- match.arg(algorithm)
  
  # check of standard tmvtnorm arguments
  cargs <- checkTmvArgs(mean, sigma, lower, upper)
  mean  <- cargs$mean
  sigma <- cargs$sigma
  lower <- cargs$lower
  upper <- cargs$upper
  
  # check of additional arguments : n and df
  if (n < 1 || !is.numeric(n) || n != as.integer(n) || length(n) > 1) {
    stop("n must be a integer scalar > 0")
  }
  
  if (df < 1 || !is.numeric(df) || length(df) > 1) {
	  stop("df must be a numeric scalar > 0")
  }
  
  if (algorithm == "rejection") {
    if (df != as.integer(df)) stop("Rejection sampling currenly works only for integer degrees of freedom. Consider using algorithm='gibbs'.")
    retval <- rtmvt.rejection(n, mean, sigma, df, lower, upper)  
  } else if (algorithm == "gibbs") {
    retval <- rtmvt.gibbs(n, mean, sigma, df, lower, upper, ...)
  }
  
  return(retval)
}

# Erzeugt eine Matrix X (n x k) mit Zufallsrealisationen aus einer Trunkierten Multivariaten t Verteilung 
# mit k Dimensionen
# ¸ber Rejection Sampling aus einer Multivariaten t-Verteilung
#
# @param n Anzahl der Realisationen
# @param mean Mittelwertvektor (k x 1) der Normalverteilung
# @param sigma Kovarianzmatrix (k x k) der Normalverteilung
# @param df degrees of freedom parameter
# @param lower unterer Trunkierungsvektor (k x 1) mit lower <= x <= upper
# @param upper oberer Trunkierungsvektor (k x 1) mit lower <= x <= upper
rtmvt.rejection <- function(n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), df = 1, 
  lower = rep(-Inf, length = length(mean)), 
  upper = rep( Inf, length = length(mean)))
{
  # No check of input parameters, checks are done in rtmvnorm()!
  
  # k = Dimension
  k <- length(mean)
  
  # mean as (1 x k) matrix
  mmean <- matrix(mean, 1, k)
  
  # Ergebnismatrix (n x k)
  Y <- matrix(NA, n, k)
  
  # Anzahl der noch zu ziehenden Samples
  numSamples <- n
  
  # Anzahl der akzeptierten Samples insgesamt
  numAcceptedSamplesTotal <- 0
  
  # Akzeptanzrate alpha aus der Multivariaten t-Verteilung bestimmen
  alpha <- pmvt(lower=lower, upper=upper, delta=mean, sigma=sigma, df=df)
  
  if (alpha <= 0.01) warning("Acceptance rate is very low and rejection sampling becomes inefficient. Consider using Gibbs sampling.")
  
  # Ziehe wiederholt aus der Multivariaten Student-t und schaue, wieviel Samples nach Trunkierung ¸brig bleiben
  while(numSamples > 0)
  {
    # Erzeuge N/alpha Samples aus einer multivariaten Normalverteilung: Wenn alpha zu niedrig ist, wird Rejection Sampling ineffizient und N/alpha zu groﬂ. Dann nur N erzeugen
    nproposals <- ifelse (numSamples/alpha > 1000000, numSamples, ceiling(max(numSamples/alpha,10)))
    X <- rmvt(nproposals, sigma=sigma, df=df) # SW: rmvt() hat keinen Parameter delta
    # add mean :  t(t(X) + mean) oder so:
    for (i in 1:k) {
      X[,i] = mean[i] + X[,i]
    }
    
    # Bestimme den Anteil der Samples nach Trunkierung
    # Bug: ind= rowSums(lower <= X & X <= upper) == k
    # wesentlich schneller als : ind=apply(X, 1, function(x) all(x >= lower & x<=upper))
    ind <- logical(nproposals)
    for (i in 1:nproposals)
    {
      ind[i] = all(X[i,] >= lower & X[i,] <= upper)
    } 
    
    # Anzahl der akzeptierten Samples in diesem Durchlauf
    numAcceptedSamples <- length(ind[ind==TRUE])
    
    # Wenn nix akzeptiert wurde, dann weitermachen
    if (length(numAcceptedSamples) == 0 || numAcceptedSamples == 0) next
    
    #cat("numSamplesAccepted=",numAcceptedSamples," numSamplesToDraw = ",numSamples,"\n")
    numNeededSamples <- min(numAcceptedSamples, numSamples)
    Y[(numAcceptedSamplesTotal+1):(numAcceptedSamplesTotal+numNeededSamples),] <- X[which(ind)[1:numNeededSamples],]
    
    # Anzahl der akzeptierten Samples insgesamt
    numAcceptedSamplesTotal <- numAcceptedSamplesTotal + numAcceptedSamples
    
    # Anzahl der verbliebenden Samples
    numSamples <- numSamples - numAcceptedSamples 
  }
  Y
}

# Gibbs sampler for the truncated multivariate Student-t
# see Geweke (1991)
#
# @param n Anzahl der Realisationen
# @param mean Mittelwertvektor (k x 1) der t-Verteilung
# @param sigma Kovarianzmatrix (k x k) der t-Verteilung
# @param df degrees of freedom parameter
# @param lower unterer Trunkierungsvektor (k x 1) mit lower <= x <= upper
# @param upper oberer Trunkierungsvektor (k x 1) mit lower <= x <= upper
# @param burn.in number of burn-in samples to be discarded
# @param start start value for Gibbs sampling
# @param thinning
rtmvt.gibbs <- function (n=1, 
  mean=rep(0, ncol(sigma)), 
  sigma = diag(length(mean)), df=1, 
  lower = rep(-Inf, length = length(mean)), 
  upper = rep( Inf, length = length(mean)),
  burn.in.samples = 0,
  start.value = NULL,
  thinning = 1)
{
  # dimension of X
  k = length(mean)    
  
  # Mean Vector
  mu = mean
  
  # number of burn-in samples
  S <- burn.in.samples
  if (!is.null(S)) {
    if (S < 0) stop("number of burn-in samples must be non-negative")   
  }
  
  # Ergebnismatrix X (n x k)
  # Random sample from truncated Student-t density 
  X <- matrix(NA, n, k)
  
  # Realisation from truncated multivariate normal
  Z <- numeric(k)
  
  # Chi-Square variable w
  w <- numeric(1)
  
  # x is one realisation from truncated Student-t density conditioned on Z and w
  x <- numeric(k)
  
  # Take start value given by user or use random start value	
  if (!is.null(start.value)) {
	  if (length(mean) != length(start.value)) stop("mean and start value have non-conforming size")
	  if (any(start.value<lower || start.value>upper)) stop("start value is not inside support region") 
	  Z <- start.value - mu 
  } else {
	  # If no start value is specified, 
    # the initial value/start value for Z drawn from TN(0,\Sigma) 
    # with truncation point a = a-mu and b = b-mu
    Z <- rtmvnorm(1, mean=rep(0,k), sigma=sigma, lower=lower-mu, upper=upper-mu, algorithm="gibbs")
  }

  # Algorithm begins :
  
  # Draw from Uni(0,1)
  U <- runif((S + n*thinning) * k)
  indU <- 1 # Index for accessing U

  # List of conditional standard deviations can be pre-calculated
  sd <- list(k)
        
  # List of t(Sigma_i) %*% solve(Sigma) term
  P  <- list(k)
  
  for(i in 1:k)
  {
    # Partitioning of Sigma
    Sigma    <- sigma[-i,-i] # (k-1) x (k-1)
    sigma_ii <- sigma[i,i]    # 1 x 1
    Sigma_i  <- sigma[i,-i]   # (k-1) x 1
    P[[i]]   <- t(Sigma_i) %*% solve(Sigma)
    sd[[i]]  <- sqrt(sigma_ii - P[[i]] %*% Sigma_i)
  }
  
  for(i in (1-S):(n*thinning))
  {
    # Step 1:   Simulation of w conditional on Z from Chi-square distribution by rejection sampling
    # so that (lower - mu) * w <= Z <= (upper - mu) * w
	  acceptedW <- FALSE
    while (!acceptedW)
    {
      w         <- (rchisq(1, df, ncp=0)/df)^(1/2)
	    acceptedW <- all((lower - mu) * w <= Z  &  Z <= (upper - mu) * w)
    }

    # Transformed Chi-Square sample subject to condition on Z0
    alpha <- (lower - mu) * w
    beta  <- (upper - mu) * w
  
    # Step 2:  Simulation from  Truncated normal Gibbs sampling approach  
    for(j in 1:k)
    {
      mu_j <- P[[j]] %*% (Z[-j])
      Fa   <- pnorm( (lower[j]-mu[j])*w, mu_j, sd[[j]])
      Fb   <- pnorm( (upper[j]-mu[j])*w, mu_j, sd[[j]])
      Z[j] <- mu_j + sd[[j]] * qnorm(U[indU] * (Fb - Fa) + Fa) #  changed on 22nd February 2010   by Manju
      indU <- indU + 1
    }
    
    #  Step 3:  Student-t transformation
    x <-  mu + ( Z / w )
    
    if (i > 0) {
      if (thinning == 1) {
        # no thinning, take all samples  except for burn-in-period
        X[i,] <- x
      }
      else if (i %% thinning == 0){
        X[i %/% thinning,] <- x  
      }
    }   
  }
  return(X)
}

# Ziehe aus einer multi-t-Distribution ohne Truncation
X <- rtmvt.rejection(n=10000, mean=rep(0, 3), df=2)

# Teste mit Kolmogoroff-Smirnoff-Test auf Verteilung
