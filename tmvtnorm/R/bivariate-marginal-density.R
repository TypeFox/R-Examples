# SW: This method is private. It is the same as mvtnorm::dmvnorm() function, 
# but without sanity checks for sigma. We perform the sanity checks before.
.dmvnorm <- function (x, mean, sigma, log = FALSE) {
  if (is.vector(x)) {
    x <- matrix(x, ncol = length(x))
  }
  distval <- mahalanobis(x, center = mean, cov = sigma)
  logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
  logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
  if (log) 
    return(logretval)
  exp(logretval)
}

# Computation of the bivariate marginal density F_{q,r}(x_q, x_r) (q != r)
# of truncated multivariate normal distribution 
# following the works of Tallis (1961), Leppard and Tallis (1989)
#
# References:
# Tallis (1961): 
#   "The Moment Generating Function of the Truncated Multi-normal Distribution"
# Leppard and Tallis (1989): 
#   "Evaluation of the Mean and Covariance of the Truncated Multinormal"
# Manjunath B G and Stefan Wilhelm (2009): 
#   "Moments Calculation for the Doubly Truncated Multivariate Normal Distribution"
#
# (n-2) Integral, d.h. zweidimensionale Randdichte in Dimension q und r, 
# da (n-2) Dimensionen rausintegriert werden.
# vgl. Tallis (1961), S.224 und Code Leppard (1989), S.550
#
# f(xq=b[q], xr=b[r])
#
# Attention: Function is not vectorized at the moment!
# Idee: Vektorisieren xq, xr --> die Integration Bounds sind immer verschieden,
#       pmvnorm() kann nicht vektorisiert werden. Sonst spart man schon ein bisschen Overhead.
# Der eigentliche bottleneck ist aber pmvnorm().
# Gibt es Unterschiede bzgl. der verschiedenen Algorithmen GenzBretz() vs. Miwa()?
# pmvnorm(corr=) kann ich verwenden
#
# @param xq
# @param xr
# @param q index for dimension q
# @param r Index für Dimension r
# @param mean
# @param sigma
# @param lower
# @param upper
# @param log=FALSE
dtmvnorm.marginal2 <- function(xq, xr, q, r, mean=rep(0, nrow(sigma)), 
  sigma=diag(length(mean)), lower=rep(-Inf, length = length(mean)), 
	upper=rep( Inf, length = length(mean)), log=FALSE,
  pmvnorm.algorithm=GenzBretz()) {
	
  # dimensionality
  n <- nrow(sigma)
  
  # number of xq values delivered
  N <- length(xq)
  
  # input checks
  if (n < 2) stop("Dimension n must be >= 2!")
  
  # TODO: Check eventuell rauslassen
  # SW; isSymmetric is sehr teuer
  #if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps))) {
  #if (!isTRUE(all.equal(sigma, t(sigma))) || any(diag(sigma) < 0)) { 
  #  stop("sigma must be a symmetric matrix")
  #}
  
  if (length(mean) != NROW(sigma)) {
    stop("mean and sigma have non-conforming size")
  }
  
  if (!(q %in% 1:n && r %in% 1:n)) {
    stop("Indexes q and r must be integers in 1:n")
  }
  
  if (q == r) {
    stop("Index q must be different than r!")
  }
  
  # Skalierungsfaktor der gestutzten Dichte (Anteil nach Trunkierung)
  # Idee: dtmvnorm.marginal2() braucht 80% der Zeit von mtmvnorm(). Die meiste Zeit davon in pmvnorm().
  # pmvnorm()-Aufrufe sind teuer, daher könnte man das alpha schon vorher berechnen
  # lassen (nur 2 pmvnorm()-Aufrufe in der Methode, würde 50% sparen)
  # Da Methode jetzt vektorisiert ist, sparen wir die Aufrufe wg. alpha
  alpha <- pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma, algorithm=pmvnorm.algorithm)
    
  if (n == 2) {
    density  <- numeric(N)
    indOut <- xq < lower[q] | xq > upper[q] | xr < lower[r] | xr > upper[r] | is.infinite(xq) | is.infinite(xr)
    density[indOut] <- 0
    # dmvnorm() macht auch viele Checks; Definiere eine private Methode .dmvnorm() ohne Checks
    density[!indOut] <- .dmvnorm(x=cbind(xq, xr)[!indOut,], mean=mean[c(q,r)], sigma=sigma[c(q,r),c(q,r)]) / alpha
    if (log == TRUE) {
	    return(log(density))
	  } else {
	    return(density)
	  }
  }
  
  # standard deviation for normalisation
  SD <- sqrt(diag(sigma))
  
  # normalised bounds
  lower.normalised  <- (lower - mean) / SD
  upper.normalised  <- (upper - mean) / SD
  
  xq.normalised     <- (xq - mean[q]) / SD[q]      # (N x 1)
  xr.normalised     <- (xr - mean[r]) / SD[r]      # (N x 1)
  
  # Computing correlation matrix R from sigma (matrix (n x n)): 
  # R = D % sigma %*% D with diagonal matrix D as sqrt(sigma)
  # same as cov2cor()
  D  <- matrix(0, n, n)
  diag(D) <- sqrt(diag(sigma))^(-1)
  R <- D %*% sigma %*% D
  
  #
  # Determine (n-2) x (n-2) correlation matrix RQR
  #
  RQR <- matrix(NA, n-2, n-2)
  RINV <- solve(R)
  WW <- matrix(NA, n-2, n-2)
  M1 <- 0
  for (i in 1:n) {
    if (i != q && i != r) {
      M1 <- M1 + 1
      M2 <- 0
      for (j in 1:n) {
        if (j != q && j != r) {
          M2 <- M2 + 1
          WW[M1, M2] <- RINV[i,j]
        }
      }
    }
  }
  WW <- solve(WW[1:(n-2),1:(n-2)])
  for(i in 1:(n-2)) {
    for(j in 1:(n-2)) {
       RQR[i, j] <- WW[i, j] / sqrt(WW[i,i] * WW[j,j])
    }
  }
  
  #
  # Determine bounds of integration vector AQR and BQR (n - 2) x 1
  #
  # lower and upper integration bounds
  AQR <- matrix(NA, N, n-2)                    
  BQR <- matrix(NA, N, n-2)
  M2 <- 0  # counter = 1..(n-2)
  for (i in 1:n) {
    if (i != q && i != r) {
      M2 <- M2 + 1
      BSQR <- (R[q, i] - R[q, r] * R[r, i]) / (1 - R[q, r]^2)    
      BSRQ <- (R[r, i] - R[q, r] * R[q, i]) / (1 - R[q, r]^2)    
      RSRQ <- (1 - R[i, q]^2) * (1 - R[q, r]^2)
      RSRQ <- (R[i, r] - R[i, q] * R[q, r]) / sqrt(RSRQ)         # partial correlation coefficient R[r,i] given q
      
      # lower integration bound
      AQR[,M2] <- (lower.normalised[i] - BSQR * xq.normalised - BSRQ * xr.normalised) / sqrt((1 - R[i, q]^2) * (1 - RSRQ^2))
      AQR[,M2] <- ifelse(is.nan(AQR[,M2]), -Inf, AQR[,M2])
      
      # upper integration bound
      BQR[,M2] <- (upper.normalised[i] - BSQR * xq.normalised - BSRQ * xr.normalised) / sqrt((1 - R[i, q]^2) * (1 - RSRQ^2))
      BQR[,M2] <- ifelse(is.nan(BQR[,M2]), Inf, BQR[,M2])
    }
  }
  
  # Correlation matrix for r and q
  R2 <- matrix(c(    1,      R[q,r], 
                R[q,r],          1), 2, 2)
              
  sigma2 <- sigma[c(q,r),c(q,r)]            
              
  density  <- ifelse (
      xq < lower[q] | 
      xq > upper[q] | 
      xr < lower[r] | 
      xr > upper[r] | is.infinite(xq) | is.infinite(xr),
     0, 
     {
      # SW: RQR is a correlation matrix, so call pmvnorm(...,corr=) which is faster than
      # pmvnorm(...,corr=)
      # SW: Possibly vectorize this loop if pmvnorm allows vectorized lower and upper bounds
      prob <- numeric(N)   # (N x 1)
      for (i in 1:N) {
        if ((n - 2) == 1) {
         # univariate case:  pmvnorm(...,corr=) does not work, will work with sigma=
         prob[i] <- pmvnorm(lower=AQR[i,], upper=BQR[i,], sigma=RQR, algorithm=pmvnorm.algorithm)
        } else {
         prob[i] <- pmvnorm(lower=AQR[i,], upper=BQR[i,], corr=RQR, algorithm=pmvnorm.algorithm)
        }
      }
	    dmvnorm(x=cbind(xq, xr), mean=mean[c(q,r)], sigma=sigma2) * prob / alpha
     }
  )
  if (log == TRUE) {
	  return(log(density))
  } else {
	  return(density)
  }
}


