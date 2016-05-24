# Estimation of the parameters 
# of the truncated multivariate normal distribution using GMM 
# and 
# (1) the moment equations from Lee (1981) and Lee (1983)
# (2) Our moment formula and equating mean and covariance matrix

#library(gmm)
#library(tmvtnorm)
#source("rtmvnorm.R")           # für checkTmvArgs()
#source("tmvnorm-estimation.R") # für vec(), vech() und inv_vech()

"%w/o%" <- function(x,y) x[!x %in% y] #--  x without y

################################################################################
#
# Multivariater Fall
#
################################################################################

# Definition einer Funktion mit Momentenbedingungen für gmm() 
# nach den Lee (1979, 1983, 1981) moment conditions
#
# N dimensions, K = N + N*(N+1)/2 parameters
# number of moment conditions L=(l_max + 1) * N
# parameter vector tet = c(mu, vech(sigma)), length K
# @param tet named parameter vector theta = c(mu, vech(sigma))
# @param x data matrix (T x N)
gmultiLee <- function(tet, fixed=c(), fullcoefnames, x, lower, upper, l_max = ceiling((ncol(x)+1)/2), cholesky=FALSE) {

 fullcoef        <- rep(NA, length(tet) + length(fixed))
 names(fullcoef) <- fullcoefnames
 if (any(!names(fixed) %in% names(fullcoef))) 
   stop("some named arguments in 'fixed' are not arguments in parameter vector theta")
 fullcoef[names(tet)]   <- tet
 fullcoef[names(fixed)] <- fixed
 
 K     <- length(tet)      # Anzahl der zu schätzenden Parameter
 N     <- ncol(x)          # Anzahl der Dimensionen
 T     <- nrow(x)          # Anzahl der Beobachtungen
 #l_max <- ceiling((N+1)/2)    # maximales l für Momentenbedingungen
 
 X     <- matrix(NA, T, (l_max+1)*N) # Rückgabematrix mit den Momenten

 # Parameter mean/sigma aus dem Parametervektor tet extrahieren
 mean <- fullcoef[1:N]
 # Matrix für sigma bauen
 if (cholesky) {
   L <- inv_vech(fullcoef[-(1:N)])
   L[lower.tri(L, diag=FALSE)] <- 0  # L entspricht jetzt chol(sigma), obere Dreiecksmatrix
   sigma <- t(L) %*% L
 } else {
   sigma <- inv_vech(fullcoef[-(1:N)])
 }
 
 #cat("Call to gmultiLee with tet=",tet," sigma=",sigma," det(sigma)=",det(sigma),"\n")
 #flush.console()
 
 # if sigma is not positive definite we return some maximum value
 if (det(sigma) <= 0 || any(diag(sigma) < 0)) {
   X     <- matrix(+Inf, T, N + N * (N+1) / 2)
   return(X)
 }
  
 sigma_inv <- solve(sigma) # inverse Kovarianzmatrix
  
 F_a     = numeric(N)
 F_b     = numeric(N)
 F       <- 1
 
 for (i in 1:N)
 {
    # one-dimensional marginal density in dimension i
    F_a[i]  <- dtmvnorm.marginal(lower[i], n=i, mean=mean, sigma=sigma, lower=lower, upper=upper)
    F_b[i]  <- dtmvnorm.marginal(upper[i], n=i, mean=mean, sigma=sigma, lower=lower, upper=upper)
 }
 
 k <- 1   
 for(l in 0:l_max) {
   for (i in 1:N)
   {
     sigma_i <- sigma_inv[i,] # i-te Zeile der inversen Kovarianzmatrix (1 x N) = entpricht sigma^{i'}
    
     a_il    <- ifelse(is.infinite(lower[i]), 0, lower[i]^l)
     b_il    <- ifelse(is.infinite(upper[i]), 0, upper[i]^l)
     
     # Lee (1983) moment equation for l
     #X[,k]   <- sigma_i %*% mean   * x[,i]^l - (x[,i]^l * x) %*% sigma_inv[,i]   + l * (x[,i]^(l-1)) + (a_il * F_a[i] - b_il * F_b[i]) / F
     X[,k]   <- sigma_i %*% mean   * x[,i]^l - sweep(x, 1, x[,i]^l, FUN="*") %*% sigma_inv[,i]   + l * (x[,i]^(l-1)) + (a_il * F_a[i] - b_il * F_b[i]) / F
     
     #T x 1     (1 x N)    (N x 1)   (T x 1)   (T x N)         (N x 1)              (T x 1)        (skalar)
     k <- k + 1 # Zählvariable
   } 
 }
 return(X)
}

# Definition einer Funktion mit Momentenbedingungen 
# mit Mean and Covariance-Matrix bauen anstatt mit Lee Bedingungen
#
# @param tet named parameter vector theta, should be part of c(vec(mu), vech(sigma))
# @param fixed a named list of fixed parameters 
# @param fullcoefnames
# @param x data matrix (T x N)
# @param lower
# @param upper
# @param cholesky flag whether we use Cholesky decompostion Sigma = LL'
#        of the covariance matrix in order to ensure positive-definiteness of sigma
gmultiManjunathWilhelm <- function(tet, fixed=c(), fullcoefnames, x, lower, upper, cholesky=FALSE) {
 
 fullcoef        <- rep(NA, length(tet) + length(fixed))
 names(fullcoef) <- fullcoefnames
 if (any(!names(fixed) %in% names(fullcoef))) 
   stop("some named arguments in 'fixed' are not arguments in parameter vector theta")
 fullcoef[names(tet)]   <- tet
 fullcoef[names(fixed)] <- fixed
 
 N     <- ncol(x)          # Anzahl der Dimensionen
 T     <- nrow(x)          # Anzahl der Beobachtungen
 
 X     <- matrix(NA, T, N + N * (N+1) / 2) # Rückgabematrix mit den Momenten

 # Parameter mean/sigma aus dem Parametervektor tet extrahieren
 mean <- fullcoef[1:N]
 # Matrix für sigma bauen
 if (cholesky) {
   L <- inv_vech(fullcoef[-(1:N)])
   L[lower.tri(L, diag=FALSE)] <- 0  # L entspricht jetzt chol(sigma), obere Dreiecksmatrix
   sigma <- t(L) %*% L
 } else {
   sigma <- inv_vech(fullcoef[-(1:N)])
 }
 
 #cat("Call to gmultiManjunathWilhelm with tet=",tet," fullcoef=", fullcoef, " sigma=",sigma," det(sigma)=",det(sigma),"\n")
 #flush.console()
 
 # if sigma is not positive definite we return some maximum value
 if (det(sigma) <= 0 || any(diag(sigma) < 0)) {
   X     <- matrix(+Inf, T, N + N * (N+1) / 2)
   return(X)
 }
 
 # Determine moments (mu, sigma) for parameters mean/sigma 
 # experimental: moments <- mtmvnorm(mean=mean, sigma=sigma, lower=lower, upper=upper, doCheckInputs=FALSE)
 moments <- mtmvnorm(mean=mean, sigma=sigma, lower=lower, upper=upper)
 
 # Momentenbedingungen für die Elemente von mean : mean(x)
 for(i in 1:N) {
   X[,i]   <- (moments$tmean[i] - x[,i])
 }
 
 # Momentenbedingungen für alle Lower-Diagonal-Elemente von sigma
 k <- 1
 for (i in 1:N) {
   for (j in 1:N) {
     # (1,1), (2, 1), (2,2)
     if (j > i) next
     #cat(sprintf("sigma[%d,%d]",i, j),"\n")
     X[,(N+k)] <- (moments$tmean[i] - x[,i]) * (moments$tmean[j] - x[,j])  - moments$tvar[i, j]
     k <- k + 1
   }
 }
 
 return(X)
}

# GMM estimation method
#
# @param X data matrix (T x n)
# @param lower, upper truncation points
# @param start list of start values for mu and sigma
# @param fixed a list of fixed parameters
# @param method either "ManjunathWilhelm" or "Lee" moment conditions
# @param cholesky flag, if TRUE, we use the Cholesky decomposition of sigma as parametrization
# @param ... additional parameters passed to gmm()
gmm.tmvnorm <- function(X, 
 lower=rep(-Inf, length = ncol(X)), 
 upper=rep(+Inf, length = ncol(X)), 
 start=list(mu=rep(0,ncol(X)), sigma=diag(ncol(X))),
 fixed=list(),
 method=c("ManjunathWilhelm","Lee"),
 cholesky=FALSE,
 ...
 ) {
 
  method <- match.arg(method)
 
  # check of standard tmvtnorm arguments
  cargs       <- checkTmvArgs(start$mu, start$sigma, lower, upper)
  start$mu    <- cargs$mean
  start$sigma <- cargs$sigma
  lower       <- cargs$lower
  upper       <- cargs$upper
  
  # check if we have at least one sample
  if (!is.matrix(X) || nrow(X) == 0) {
    stop("Data matrix X with at least one row required.")
  }
  
  # verify dimensions of x and lower/upper match
  n <- length(lower)
  if (NCOL(X) != n) {
		stop("data matrix X has a non-conforming size. Must have ",length(lower)," columns.")
	}
	
	# check if lower <= X <= upper for all rows
	ind <- logical(nrow(X))
  for (i in 1:nrow(X))
  {
    ind[i] = all(X[i,] >= lower & X[i,] <= upper)
  }
  if (!all(ind)) {
    stop("some of the data points are not in the region lower <= X <= upper")
  }
  
  # parameter vector theta
  theta <- c(start$mu, vech2(start$sigma)) 
  # names for mean vector elements : mu_i
  nmmu     <- paste("mu_",1:n,sep="")
  # names for sigma elements : sigma_i.j
  nmsigma  <- paste("sigma_",vech2(outer(1:n,1:n, paste, sep=".")),sep="")
  names(theta)  <- c(nmmu, nmsigma)
  fullcoefnames <- names(theta)
  
  # use only those parameters without the fixed parameters for gmm(),
  # since I do not know how to specify fixed=c() in gmm()
  theta2 <- theta[names(theta) %w/o% names(fixed)]
  
  # define a wrapper function with only 2 arguments theta and x (f(theta, x)) 
  # that will be invoked by gmm()
  gManjunathWilhelm <- function(tet, x) {
    gmultiManjunathWilhelm(tet=tet, fixed=unlist(fixed), 
      fullcoefnames=fullcoefnames, x=x, 
      lower=lower, upper=upper, cholesky=cholesky)
  }
  
  # TODO: Allow for l_max parameter for Lee moment conditions
  gLee <- function(tet, x) {
    gmultiLee(tet = tet, fixed = unlist(fixed), 
      fullcoefnames = fullcoefnames, x = x, 
      lower = lower, upper = upper, cholesky = cholesky)
  }
  
  if (method == "ManjunathWilhelm") {
    gmm.fit <- gmm(gManjunathWilhelm, x=X, t0=theta2, ...)
  } else {
    gmm.fit <- gmm(gLee, x=X, t0=theta2, ...)
  }
  return(gmm.fit)
}

# deprecated
# GMM mit Lee conditions
gmm.tmvnorm2 <- function (X, lower = rep(-Inf, length = ncol(X)), upper = rep(+Inf, 
    length = ncol(X)), start = list(mu = rep(0, ncol(X)), sigma = diag(ncol(X))), 
    fixed = list(), cholesky = FALSE, ...) 
{
    cargs <- checkTmvArgs(start$mu, start$sigma, lower, upper)
    start$mu <- cargs$mean
    start$sigma <- cargs$sigma
    lower <- cargs$lower
    upper <- cargs$upper
    if (!is.matrix(X) || nrow(X) == 0) {
        stop("Data matrix X with at least one row required.")
    }
    n <- length(lower)
    if (NCOL(X) != n) {
        stop("data matrix X has a non-conforming size. Must have ", 
            length(lower), " columns.")
    }
    ind <- logical(nrow(X))
    for (i in 1:nrow(X)) {
        ind[i] = all(X[i, ] >= lower & X[i, ] <= upper)
    }
    if (!all(ind)) {
        stop("some of the data points are not in the region lower <= X <= upper")
    }
    theta <- c(start$mu, vech2(start$sigma))
    nmmu <- paste("mu_", 1:n, sep = "")
    nmsigma <- paste("sigma_", vech2(outer(1:n, 1:n, paste, sep = ".")), 
        sep = "")
    names(theta) <- c(nmmu, nmsigma)
    fullcoefnames <- names(theta)
    theta2 <- theta[names(theta) %w/o% names(fixed)]
    gmultiwrapper <- function(tet, x) {
        gmultiLee(tet = tet, fixed = unlist(fixed), 
            fullcoefnames = fullcoefnames, x = x, lower = lower, 
            upper = upper, cholesky = cholesky)
    }
    gmm.fit <- gmm(gmultiwrapper, x = X, t0 = theta2, ...)
    return(gmm.fit)
}




