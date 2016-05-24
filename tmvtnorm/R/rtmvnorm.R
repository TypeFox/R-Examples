################################################################################
#
# Sampling from Truncated multivariate Gaussian distribution using 
#
# a) Rejection sampling
# b) Gibbs sampler
#
# for both rectangular constraints a <= x <= b and general linear constraints
# a <= Dx <= b. For D = I this implies rectangular constraints.
# The method can be used using both covariance matrix sigma and precision matrix H.
# 
# Author: Stefan Wilhelm
#
# References:
# (1) Jayesh H. Kotecha and Petar M. Djuric (1999) : 
# "GIBBS SAMPLING APPROACH FOR GENERATION OF TRUNCATED MULTIVARIATE GAUSSIAN RANDOM VARIABLES"
# (2) Geweke (1991): 
# "Effcient simulation from the multivariate normal and Student-t distributions 
# subject to linear constraints and the evaluation of constraint probabilities"
# (3) John Geweke (2005): Contemporary Bayesian Econometrics and Statistics, Wiley, pp.171-172
# (4) Wilhelm (2011) package vignette to package "tmvtnorm"
#
################################################################################

# We need this separate method rtmvnorm.sparseMatrix() because
# rtmvnorm() initialises dense d x d sigma and D matrix which will not work for high dimensions d.
# It also does some sanity checks on sigma and D (determinant etc.) which will not
# work for high dimensions.
 
# returns a matrix X (n x d) with random draws 
# from a truncated multivariate normal distribution with d dimensionens
# using Gibbs sampling
#
# @param n Anzahl der Realisationen
# @param mean mean vector (d x 1) der Normalverteilung
# @param lower lower truncation vector (d x 1) with lower <= x <= upper
# @param upper upper truncation vector (d x 1) with lower <= x <= upper
# @param H precision matrix (d x d) if given, defaults to identity matrix
rtmvnorm.sparseMatrix <- function(n, 
    mean = rep(0, nrow(H)), 
    H = sparseMatrix(i=1:length(mean), j=1:length(mean), x=1),
    lower = rep(-Inf, length = length(mean)), 
    upper = rep( Inf, length = length(mean)),
    ...)
{
  if (is.null(H) || !inherits(H, "sparseMatrix")) {
    stop("H must be of class 'sparseMatrix'")
  }
  rtmvnorm.gibbs.Fortran(n, mean, sigma=NULL, H, lower, upper, ...)
}

# Erzeugt eine Matrix X (n x d) mit Zufallsrealisationen 
# aus einer Trunkierten Multivariaten Normalverteilung mit d Dimensionen
# ¸ber Rejection Sampling oder Gibbs Sampler aus einer Multivariaten Normalverteilung.
# If matrix D is given, it must be a (d x d) full rank matrix.
# Therefore this method can only cover the case with only r <= d linear restrictions.
# For r > d linear restrictions, please see rtmvnorm2(n, mean, sigma, D, lower, upper),
# where D can be defined as (r x d).
#
# @param n Anzahl der Realisationen
# @param mean Mittelwertvektor (d x 1) der Normalverteilung
# @param sigma Kovarianzmatrix (d x d) der Normalverteilung
# @param lower unterer Trunkierungsvektor (d x 1) mit lower <= Dx <= upper
# @param upper oberer Trunkierungsvektor (d x 1) mit lower <= Dx <= upper
# @param D Matrix for linear constraints, defaults to (d x d) diagonal matrix
# @param H Precision matrix (d x d) if given
# @param algorithm c("rejection", "gibbs", "gibbsR")
rtmvnorm <- function(n, 
    mean = rep(0, nrow(sigma)), 
    sigma = diag(length(mean)),
    lower = rep(-Inf, length = length(mean)), 
    upper = rep( Inf, length = length(mean)),
    D = diag(length(mean)),
    H = NULL, 
    algorithm=c("rejection", "gibbs", "gibbsR"), ...)
{
  algorithm <- match.arg(algorithm)
  
  if (is.null(mean) && (is.null(sigma) || is.null(H))) {
    stop("Invalid arguments for ",sQuote("mean")," and ",sQuote("sigma"),"/",sQuote("H"),". Need at least mean vector and covariance or precision matrix.")
  }
  
  # check of standard tmvtnorm arguments
  cargs <- checkTmvArgs(mean, sigma, lower, upper)
  mean  <- cargs$mean
  sigma <- cargs$sigma
  lower <- cargs$lower
  upper <- cargs$upper
  
  if (!is.null(H) && sigma != diag(length(mean))) {
   stop("Cannot give both covariance matrix sigma and precision matrix H arguments at the same time")
  }
  else if (!is.null(H) && !inherits(H, "sparseMatrix")) {
   # check precision matrix H if it is symmetric and positive definite
   checkSymmetricPositiveDefinite(H, name="H")
   # H explicitly given, we will override sigma later if we need sigma
   # sigma <- solve(H)
  } 
  # else sigma explicitly or implicitly given
    
  # check of additional arguments
  if (n < 1 || !is.numeric(n) || n != as.integer(n) || length(n) > 1) {
	  stop("n must be a integer scalar > 0")
  }
  
  # check matrix D, must be n x n with rank n
  if (!is.matrix(D) || det(D) == 0) {
    stop("D must be a (n x n) matrix with full rank n!")
  }
  
  if (!identical(D,diag(length(mean)))) {
    # D <> I : general linear constraints
    retval <- rtmvnorm.linear.constraints(n=n, mean=mean, sigma=sigma, H=H, lower=lower, upper=upper, D=D, algorithm=algorithm, ...)
    return(retval)
  } else {
    # D == I : rectangular case
    if (algorithm == "rejection") {
      if (!is.null(H)) {
        # precision matrix case H
        retval <- rtmvnorm.rejection(n, mean, sigma=solve(H), lower, upper, ...)
      } else {
        # covariance matrix case sigma
        retval <- rtmvnorm.rejection(n, mean, sigma, lower, upper, ...)
      }
    } else if (algorithm == "gibbs") {
      # precision matrix case H vs. covariance matrix case sigma will be handled inside method 
      retval <- rtmvnorm.gibbs.Fortran(n, mean, sigma, H, lower, upper, ...)
    } else if (algorithm == "gibbsR") {
      if (!is.null(H)) {
        # precision matrix case H
        retval <- rtmvnorm.gibbs.Precision(n, mean, H, lower, upper, ...)  
      } else {
        # covariance matrix case sigma
        retval <- rtmvnorm.gibbs(n, mean, sigma, lower, upper, ...)  
      }
    }
  }
  return(retval)
}

# Erzeugt eine Matrix X (n x k) mit Zufallsrealisationen 
# aus einer Trunkierten Multivariaten Normalverteilung mit k Dimensionen
# ¸ber Rejection Sampling aus einer Multivariaten Normalverteilung mit der Bedingung
# lower <= Dx <= upper
# 
# Wenn D keine Diagonalmatrix ist, dann ist gelten lineare Restriktionen f¸r
# lower <= Dx <= upper (siehe Geweke (1991))
#
# @param n Anzahl der Realisationen
# @param mean Mittelwertvektor (k x 1) der Normalverteilung
# @param sigma Kovarianzmatrix (k x k) der Normalverteilung
# @param lower unterer Trunkierungsvektor (k x 1) mit lower <= x <= upper
# @param upper oberer Trunkierungsvektor (k x 1) mit lower <= x <= upper
# @param D Matrix for linear constraints, defaults to diagonal matrix
rtmvnorm.rejection <- function(n, 
  mean = rep(0, nrow(sigma)), 
  sigma = diag(length(mean)), 
  lower = rep(-Inf, length = length(mean)), 
  upper = rep( Inf, length = length(mean)),
  D = diag(length(mean)))
{
  # No check of input parameters, checks are done in rtmvnorm()!
  
  # k = Dimension
  k <- length(mean)
  
  # Ergebnismatrix (n x k)
  Y <- matrix(NA, n, k)
  
  # Anzahl der noch zu ziehenden Samples
  numSamples <- n
  
  # Anzahl der akzeptierten Samples insgesamt
  numAcceptedSamplesTotal <- 0
  
  # Akzeptanzrate alpha aus der Multivariaten Normalverteilung bestimmen
  r <- length(lower)
  d <- length(mean)
  if (r == d & identical(D, diag(d))) {
    alpha <- pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma)
    if (alpha <= 0.01) warning(sprintf("Acceptance rate is very low (%s) and rejection sampling becomes inefficient. Consider using Gibbs sampling.", alpha))
    estimatedAlpha <- TRUE
  } else {
    # TODO: Wie bestimme ich aus lower <= Dx <= upper f¸r r > d Restriktionen die Akzeptanzrate alpha?
    # Defere calculation of alpha. Assume for now that all samples will be accepted.
    alpha <- 1
    estimatedAlpha <- FALSE
  }
  
  # Ziehe wiederholt aus der Multivariaten NV und schaue, wieviel Samples nach Trunkierung ¸brig bleiben
  while(numSamples > 0)
  {
    # Erzeuge N/alpha Samples aus einer multivariaten Normalverteilung: Wenn alpha zu niedrig ist, wird Rejection Sampling ineffizient und N/alpha zu groﬂ. Dann nur N erzeugen
    nproposals <- ifelse (numSamples/alpha > 1000000, numSamples, ceiling(max(numSamples/alpha,10)))
    X <- rmvnorm(nproposals, mean=mean, sigma=sigma)
    
    # Bestimme den Anteil der Samples nach Trunkierung
    # Bug: ind= rowSums(lower <= X & X <= upper) == k
    # wesentlich schneller als : ind=apply(X, 1, function(x) all(x >= lower & x<=upper))
    X2 <- X %*% t(D)
    ind <- logical(nproposals)
    for (i in 1:nproposals)
    {
      ind[i] <- all(X2[i,] >= lower & X2[i,] <= upper)
    } 

    # Anzahl der akzeptierten Samples in diesem Durchlauf
    numAcceptedSamples <- length(ind[ind==TRUE])
        
    # Wenn nix akzeptiert wurde, dann weitermachen
    if (length(numAcceptedSamples) == 0 || numAcceptedSamples == 0) next
    
    if (!estimatedAlpha) {
      alpha <- numAcceptedSamples / nproposals
      if (alpha <= 0.01) warning(sprintf("Acceptance rate is very low (%s) and rejection sampling becomes inefficient. Consider using Gibbs sampling.", alpha))
    }
    
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

# Gibbs Sampler for Truncated Univariate Normal Distribution
#
# Jayesh H. Kotecha and Petar M. Djuric (1999) : GIBBS SAMPLING APPROACH FOR GENERATION OF TRUNCATED MULTIVARIATE GAUSSIAN RANDOM VARIABLES
#
# Im univariaten Fall sind die erzeugten Samples unabh‰ngig, 
# deswegen gibt es hier keine Chain im eigentlichen Sinn und auch keinen Startwert/Burn-in/Thinning. 
#
# As a change to Kotecha, we do not draw a sample x from the Gaussian Distribution
# and then apply pnorm(x) - which is uniform - but rather draw directly from the
# uniform distribution u ~ U(0, 1).
#
# @param n number of realisations
# @param mu      mean of the normal distribution
# @param sigma   standard deviation
# @param a lower truncation point
# @param b upper truncation point
rtnorm.gibbs <- function(n, mu=0, sigma=1, a=-Inf, b=Inf)
{
   # Draw from Uni(0,1)
   F <- runif(n) 	
   
   #Phi(a) und Phi(b)
   Fa <- pnorm(a, mu, sd=sigma)
   Fb <- pnorm(b, mu, sd=sigma)
   
   # Truncated Normal Distribution, see equation (6), but F(x) ~ Uni(0,1), 
   # so we directly draw from Uni(0,1) instead of doing:
   # x <- rnorm(n, mu, sigma)
   # y <-  mu + sigma * qnorm(pnorm(x)*(Fb - Fa) + Fa)
   y  <-  mu + sigma * qnorm(F * (Fb - Fa) + Fa)	
   
   y
}

# Gibbs Sampler Implementation in R for Truncated Multivariate Normal Distribution
# (covariance case with sigma)
# Jayesh H. Kotecha and Petar M. Djuric (1999) : 
# GIBBS SAMPLING APPROACH FOR GENERATION OF TRUNCATED MULTIVARIATE 
# GAUSSIAN RANDOM VARIABLES
#
#
# @param n Anzahl der Realisationen
# @param mean Mittelwertvektor (k x 1) der Normalverteilung
# @param sigma Kovarianzmatrix (k x k) der Normalverteilung
# @param lower unterer Trunkierungsvektor (k x 1) mit lower <= Dx <= upper
# @param upper oberer Trunkierungsvektor (k x 1) mit lower <= Dx <= upper
# @param burn.in number of burn-in samples to be discarded
# @param start start value for Gibbs sampling
# @param thinning
rtmvnorm.gibbs <- function(n, 
    mean = rep(0, nrow(sigma)), 
    sigma = diag(length(mean)),
    lower = rep(-Inf, length = length(mean)), 
    upper = rep( Inf, length = length(mean)), 
    burn.in.samples = 0, 
    start.value = NULL, 
    thinning = 1)
{
  # We check only additional arguments like "burn.in.samples", "start.value" and "thinning"
  
  if (thinning < 1 || !is.numeric(thinning) || length(thinning) > 1) {
	  stop("thinning must be a integer scalar > 0")
  }
	
  # dimension of X
  d <- length(mean)
  
  # number of burn-in samples
  S <- burn.in.samples
  if (!is.null(S)) {
	if (S < 0) stop("number of burn-in samples must be non-negative")   
  }
	
  # Take start value given by user or determine from lower and upper	
  if (!is.null(start.value)) {
    if (length(mean) != length(start.value)) stop("mean and start value have non-conforming size")
	  if (any(start.value<lower || start.value>upper)) stop("start value is not inside support region") 
	  x0 <- start.value 
  } else {
    # Start value from support region, may be lower or upper bound, if they are finite, 
	  # if both are infinite, we take 0.
	  x0  <- ifelse(is.finite(lower), lower, ifelse(is.finite(upper), upper, 0))
  }
  
  # Sample from univariate truncated normal distribution which is very fast.
  if (d == 1)
  {
    X <- rtnorm.gibbs(n, mu=mean[1], sigma=sigma[1,1], a=lower[1], b=upper[1])
    return(X)
  }
      
  # Ergebnismatrix (n x k)
  X <- matrix(NA, n, d)
  
  # Draw from Uni(0,1)
  U <- runif((S + n*thinning) * d)
  l <- 1
  
  # List of conditional standard deviations can be pre-calculated
  sd <- list(d)
  # List of t(Sigma_i) %*% solve(Sigma) term
  P  <- list(d)
  
  for(i in 1:d)
  {
    # Partitioning of Sigma
    Sigma    <- sigma[-i,-i] # (d-1) x (d-1)
    sigma_ii <- sigma[i,i]   # 1 x 1
    Sigma_i  <- sigma[i,-i]  # 1 x (d-1)
    
    P[[i]]   <- t(Sigma_i) %*% solve(Sigma)  # (1 x (d-1)) * ((d-1) x (d-1)) =  (1 x (d-1))
    sd[[i]]  <- sqrt(sigma_ii - P[[i]] %*% Sigma_i)  # (1 x (d-1)) * ((d-1) x 1)
  }
  
  x <- x0
  
  # Runn chain from index (1 - #burn-in-samples):(n*thinning) and only record samples from j >= 1
  # which discards the burn-in-samples
  for (j in (1-S):(n*thinning))
  {
    # For all dimensions
    for(i in 1:d)
    {
      # Berechnung von bedingtem Erwartungswert und bedingter Varianz:
      # bedingte Varianz h‰ngt nicht von x[-i] ab!
      mu_i  <- mean[i]    + P[[i]] %*% (x[-i] - mean[-i])
      
      # Transformation
	    F.tmp <- pnorm(c(lower[i], upper[i]), mu_i, sd[[i]])
	    Fa    <- F.tmp[1]
      Fb    <- F.tmp[2]
	    x[i]  <- mu_i + sd[[i]] * qnorm(U[l] * (Fb - Fa) + Fa)
      l     <- l + 1
    }
	
	  if (j > 0) {
	    if (thinning == 1) {
	      # no thinning, take all samples	except for burn-in-period
	      X[j,] <- x
	    }
	    else if (j %% thinning == 0){
	      X[j %/% thinning,] <- x	
	    }
    }
  }
  return(X)
}

# R-Implementation of Gibbs sampler with precision matrix H
#
# @param n number of random draws
# @param mean Mittelwertvektor (k x 1) der Normalverteilung
# @param H Precision matrix (k x k) der Normalverteilung
# @param lower unterer Trunkierungsvektor (k x 1) mit lower <= x <= upper
# @param upper oberer Trunkierungsvektor (k x 1) mit lower <= x <= upper
# @param burn.in number of burn-in samples to be discarded
# @param start start value for Gibbs sampling
# @param thinning
rtmvnorm.gibbs.Precision <- function(n,
    mean = rep(0, nrow(H)),
    H = diag(length(mean)),
    lower = rep(-Inf, length = length(mean)),
    upper = rep( Inf, length = length(mean)),
    burn.in.samples = 0,
    start.value = NULL,
    thinning = 1)
{
  # We check only additional arguments like "burn.in.samples", "start.value" and "thinning"
  if (thinning < 1 || !is.numeric(thinning) || length(thinning) > 1) {
	  stop("thinning must be a integer scalar > 0")
  }

  # dimension of X
  d <- length(mean)

  # number of burn-in samples
  S <- burn.in.samples
  if (!is.null(S)) {
	  if (S < 0) stop("number of burn-in samples must be non-negative")
  }

  # Take start value given by user or determine from lower and upper
  if (!is.null(start.value)) {
    if (length(mean) != length(start.value)) stop("mean and start value have non-conforming size")
	  if (any(start.value<lower || start.value>upper)) stop("start value is not inside support region")
	  x0 <- start.value
  } else {
    # Start value from support region, may be lower or upper bound, if they are finite,
	  # if both are infinite, we take 0.
	  x0  <- ifelse(is.finite(lower), lower, ifelse(is.finite(upper), upper, 0))
  }

  # Sample from univariate truncated normal distribution which is very fast.
  if (d == 1) {
    X <- rtnorm.gibbs(n, mu=mean[1], sigma=1/H[1,1], a=lower[1], b=upper[1])
    return(X)
  }

  # Ergebnismatrix (n x k)
  X <- matrix(NA, n, d)

  # Draw from U ~ Uni(0,1) for all iterations we need in advance
  U <- runif((S + n*thinning) * d)
  l <- 1

  # Vector of conditional standard deviations sd(i | -i) = H_ii^{-1} = 1 / H[i, i] = sqrt(1 / diag(H))
  # does not depend on x[-i] and can be precalculated before running the chain.
  sd <- sqrt(1 / diag(H))
  
  # start value of the chain
  x <- x0

  # Run chain from index (1 - #burn-in-samples):(n*thinning) and only record samples from j >= 1
  # which discards the burn-in-samples
  for (j in (1-S):(n*thinning))
  {
    # For all dimensions
    for(i in 1:d)
    {
      # conditional mean mu[i] = E[i | -i] = mean[i] - H_ii^{-1} H[i,-i] (x[-i] - mean[-i])
      mu_i  <- mean[i]    - (1 / H[i,i]) * H[i,-i] %*% (x[-i] - mean[-i])

      # draw x[i | -i] from conditional univariate truncated normal distribution with
      # TN(E[i | -i], sd(i | -i), lower[i], upper[i]) 
	    F.tmp <- pnorm(c(lower[i], upper[i]), mu_i, sd[i])
	    Fa    <- F.tmp[1]
      Fb    <- F.tmp[2]
	    x[i]  <- mu_i + sd[i] * qnorm(U[l] * (Fb - Fa) + Fa)
      l     <- l + 1
    }

	  if (j > 0) {
	    if (thinning == 1) {
	     # no thinning, take all samples	except for burn-in-period
	     X[j,] <- x
	    }
	    else if (j %% thinning == 0){
	     X[j %/% thinning,] <- x
	    }
    }
  }
  return(X)
}

# Gibbs sampler with compiled Fortran code
# Depending on, whether covariance matrix Sigma or precision matrix H (dense or sparse format) 
# is specified as parameter, we call either 
# Fortran routine "rtmvnormgibbscov" (dense covariance matrix sigma), 
# "rtmvnormgibbsprec" (dense matrix H) or "rtmvnormgibbssparseprec" (sparse precision matrix H).
#
# @param H precision matrix in sparse triplet format (i, j, v)
# Memory issues: We want to increase dimension d, and return matrix X will be (n x d)
# so if we want to create a large number of random samples X (n x d) with high d then
# we will probably also run into memory problems (X is dense). In most MCMC applications,
# we only have to create a small number n in high dimensions, 
# e.g. 1 random sample per iteration (+ burn-in-samples). 
# In this case we will not experience any problems. Users should be aware of choosing n and d appropriately
rtmvnorm.gibbs.Fortran <- function(n, 
    mean = rep(0, nrow(sigma)), 
    sigma = diag(length(mean)),
    H     = NULL, 
    lower = rep(-Inf, length = length(mean)), 
    upper = rep( Inf, length = length(mean)), 
		burn.in.samples = 0, start.value = NULL, thinning = 1)
{
  # No checks of input arguments, checks are done in rtmvnorm()
  
  # dimension of X
  d <- length(mean)
  
  # number of burn-in samples
  S <- burn.in.samples
  if (!is.null(S)) {
	if (S < 0) stop("number of burn-in samples must be non-negative")   
  }
	
  # Take start value given by user or determine from lower and upper	
  if (!is.null(start.value)) {
    if (length(mean) != length(start.value)) stop("mean and start value have non-conforming size")
	if (any(start.value<lower || start.value>upper)) stop("start value is not inside support region") 
	x0 <- start.value 
  } else {
    # Start value from support region, may be lower or upper bound, if they are finite, 
	# if both are infinite, we take 0.
	x0  <- ifelse(is.finite(lower), lower, ifelse(is.finite(upper), upper, 0))
  }
  
  # Sample from univariate truncated normal distribution which is very fast.
  if (d == 1) {
    if (!is.null(H)) {
      X <- rtnorm.gibbs(n, mu=mean[1], sigma=1 / sigma[1,1], a=lower[1], b=upper[1])
    } else {
      X <- rtnorm.gibbs(n, mu=mean[1], sigma=sigma[1,1], a=lower[1], b=upper[1])
    }
    return(X)
  }
      
  # Ergebnismatrix (n x d)
  X <- matrix(0, n, d)
  
  # Call to Fortran subroutine
  if (!is.null(H)){
    if (!inherits(H, "sparseMatrix")) {
      ret <- .Fortran("rtmvnormgibbsprec",
                              n     = as.integer(n),
                              d     = as.integer(d),
                              mean  = as.double(mean),
                              H     = as.double(H),
                              lower = as.double(lower), 
                              upper = as.double(upper),
                              x0    = as.double(x0),
							  burnin   = as.integer(burn.in.samples),
							  thinning = as.integer(thinning),
                              X     = as.double(X), 
                              NAOK=TRUE, PACKAGE="tmvtnorm")
    } else if (inherits(H, "dgCMatrix")) {  # H is given in compressed sparse column (csc) representation
      ret <- .Fortran("rtmvnorm_sparse_csc",
                              n     = as.integer(n),
                              d     = as.integer(d),
                              mean  = as.double(mean),
                              Hi    = as.integer(H@i),
                              Hp    = as.integer(H@p),
                              Hv    = as.double(H@x),
                              num_nonzero = as.integer(length(H@x)),
                              lower = as.double(lower), 
                              upper = as.double(upper),
                              x0    = as.double(x0),
							                burnin   = as.integer(burn.in.samples),
							                thinning = as.integer(thinning),
                              X     = as.double(X), 
                              NAOK=TRUE, PACKAGE="tmvtnorm")
    }
    else { # H is given in sparse matrix triplet representation
      # Es muss klar sein, dass nur die obere Dreiecksmatrix (i <= j) ¸bergeben wird...
      sH <- as(H, "dgTMatrix")        # precision matrix as triplet representation 
      # ATTENTION: sH@i and sH@j are zero-based (0..(n-1)), we need it as 1...n
      ind <- sH@i <= sH@j             # upper triangular matrix elements of H[i,j] with i <= j
      ret <- .Fortran("rtmvnorm_sparse_triplet",
                              n     = as.integer(n),
                              d     = as.integer(d),
                              mean  = as.double(mean),
                              Hi    = as.integer(sH@i[ind]+1),
                              Hj    = as.integer(sH@j[ind]+1),
                              Hv    = as.double(sH@x[ind]),
                              num_nonzero = as.integer(sum(ind)),
                              lower = as.double(lower), 
                              upper = as.double(upper),
                              x0    = as.double(x0),
							                burnin   = as.integer(burn.in.samples),
							                thinning = as.integer(thinning),
                              X     = as.double(X), 
                              NAOK=TRUE, PACKAGE="tmvtnorm")
    }
  } else {
    ret <- .Fortran("rtmvnormgibbscov",
                              n     = as.integer(n),
                              d     = as.integer(d),
                              mean  = as.double(mean),
                              sigma = as.double(sigma),
                              lower = as.double(lower), 
                              upper = as.double(upper),
                              x0    = as.double(x0),
							  burnin   = as.integer(burn.in.samples),
							  thinning = as.integer(thinning),
                              X     = as.double(X), 
                              NAOK=TRUE, PACKAGE="tmvtnorm")
  }
  X <- matrix(ret$X, ncol=d, byrow=TRUE)
  return(X)
}


# Gibbs sampling f¸r Truncated Multivariate Normal Distribution 
# with linear constraints based on Geweke (1991): 
# This is simply a wrapper function around our rectangular sampling version...
#
# x ~ N(mu, sigma) subject to a <= Dx <= b
#
# alpha <= z <= beta 
# mit alpha = a - D * mu, beta = b - D * mu
# z ~ N(0, T), T = D Sigma D'
# x = mu + D^(-1) z
#
# @param n Anzahl der Realisationen
# @param mean Mittelwertvektor (k x 1) der t-verteilung
# @param sigma Kovarianzmatrix (k x k) der t-Verteilung
# @param lower unterer Trunkierungsvektor (k x 1) mit lower <= x <= upper
# @param upper oberer Trunkierungsvektor (k x 1) mit lower <= x <= upper
# @param D Matrix for linear constraints, defaults to diagonal matrix
# @param burn.in number of burn-in samples to be discarded
# @param start start value for Gibbs sampling
# @param thinning
rtmvnorm.linear.constraints <- 
  function(n, 
  mean = rep(0, nrow(sigma)), 
  sigma = diag(length(mean)),
  H = NULL, 
  lower = rep(-Inf, length = length(mean)), 
  upper = rep( Inf, length = length(mean)),
  D  = diag(length(mean)), 
  algorithm,...)
{
  # dimension of X
  d <- length(mean)
  
  # check matrix D, must be n x n with rank n
  if (!is.matrix(D) || det(D) == 0) {
    stop("D must be a (n x n) matrix with full rank n!")
  }
  
  # create truncated multi-normal samples in variable Z ~ N(0, T) 
  # with alpha <= z <= beta 
	
  # Parameter-Transformation for given sigma:
  # x ~ N(mean, sigma) subject to a <= Dx <= b
  # define z = D x - D mu
  # alpha <= z <= beta 
  # mit alpha = a - D * mu
  #     beta  = b - D * mu
  # z ~ N(0, T), 
  # T = D Sigma D'
  # x = mu + D^(-1) z
  # Parameter-Transformation for given H:
  # x ~ N(mean, H^{-1})
  # precision matrix in z is:
  # T^{-1} = D'^{-1} H D^{-1}    # (AB)^{-1} = B^{-1} %*% A^{-1}
  alpha <- as.vector(lower - D %*% mean)
  beta  <- as.vector(upper - D %*% mean)
  Dinv  <- solve(D) # D^(-1)
  
  if (!is.null(H)) {
    Tinv <- t(Dinv) %*% H %*% Dinv
    Z <- rtmvnorm(n, mean=rep(0, d), sigma=diag(d), H=Tinv, lower=alpha, upper=beta, algorithm=algorithm, ...)
  } else {
    T     <- D %*% sigma %*% t(D)
	  Z <- rtmvnorm(n, mean=rep(0, d), sigma=T, H=NULL, lower=alpha, upper=beta, algorithm=algorithm, ...)
  }

  # For each z do the transformation
  # x = mu + D^(-1) z
  X <- sweep(Z %*% t(Dinv), 2, FUN="+", mean)
  return(X)
}

################################################################################

if (FALSE) {

checkSymmetricPositiveDefinite(matrix(1:4, 2, 2), name = "H")

lower <- c(-1, -1)
upper <- c(1, 1)
mean <- c(0.5, 0.5)
sigma <- matrix(c(1, 0.8, 0.8, 1), 2, 2)
H <- solve(sigma)
D <- matrix(c(1, 1, 1, -1), 2, 2)

checkSymmetricPositiveDefinite(H, name = "H")

# 1. covariance matrix sigma case
# 1.1. rectangular case D == I
X0 <- rtmvnorm(n=1000, mean, sigma, lower, upper, algorithm="rejection")
X1 <- rtmvnorm(n=1000, mean=mean, sigma=sigma, lower=lower, upper=upper, algorithm="rejection")
X2 <- rtmvnorm(n=1000, mean=mean, sigma=sigma, lower=lower, upper=upper, algorithm="gibbsR")
X3 <- rtmvnorm(n=1000, mean=mean, sigma=sigma, lower=lower, upper=upper, algorithm="gibbs")

par(mfrow=c(2,2))
plot(X1)
plot(X2)
plot(X3)

cov(X1)
cov(X2)
cov(X3)

# 1.2. general linear constraints case D <> I
X1 <- rtmvnorm(n=1000, mean=mean, sigma=sigma, lower=lower, upper=upper, D=D, algorithm="rejection")
X2 <- rtmvnorm(n=1000, mean=mean, sigma=sigma, lower=lower, upper=upper, D=D, algorithm="gibbsR")
X3 <- rtmvnorm(n=1000, mean=mean, sigma=sigma, lower=lower, upper=upper, D=D, algorithm="gibbs")

par(mfrow=c(2,2))
plot(X1)
plot(X2)
plot(X3)

# 2. precision matrix case H
# 2.1. rectangular case D == I
X1 <- rtmvnorm(n=1000, mean=mean, H=H, lower=lower, upper=upper, algorithm="rejection")
X2 <- rtmvnorm(n=1000, mean=mean, H=H, lower=lower, upper=upper, algorithm="gibbsR")
X3 <- rtmvnorm(n=1000, mean=mean, H=H, lower=lower, upper=upper, algorithm="gibbs")

par(mfrow=c(2,2))
plot(X1)
plot(X2)
plot(X3)

# 2.2. general linear constraints case D <> I
X1 <- rtmvnorm(n=1000, mean=mean, H=H, lower=lower, upper=upper, D=D, algorithm="rejection")
X2 <- rtmvnorm(n=1000, mean=mean, H=H, lower=lower, upper=upper, D=D, algorithm="gibbsR")
X3 <- rtmvnorm(n=1000, mean=mean, H=H, lower=lower, upper=upper, D=D, algorithm="gibbs")

par(mfrow=c(2,2))
plot(X1)
plot(X2)
plot(X3)

}
