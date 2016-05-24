# Expectation and covariance matrix computation 
# based on the algorithms by Lee (1979), Lee (1983), Leppard and Tallis (1989)
# and Manjunath and Wilhelm (2009)
#
# References:
# Amemiya (1973) : Regression Analysis When the Dependent Variable is Truncated Normal
# Amemiya (1974) : Multivariate Regression and Simultaneous Equations Models When the Dependent Variables Are Truncated Normal
# Lee (1979)     : On the first and second moments of the truncated multi-normal distribution and a simple estimator
# Lee (1983)     : The Determination of Moments of the Doubly Truncated Multivariate Tobit Model
# Leppard and Tallis (1989) : Evaluation of the Mean and Covariance of the Truncated Multinormal
# Manjunath B G and Stefan Wilhelm (2009): 
#   Moments Calculation for the Doubly Truncated Multivariate Normal Distribution
# Johnson/Kotz (1972)

# Compute truncated mean and truncated variance in the case 
# where only a subset of k < n  x_1,..,x_k variables are truncated.
# In this case, computations simplify and we only have to deal with k-dimensions.
# Example: n=10 variables but only k=3 variables are truncated.
#
# Attention: Johnson/Kotz (1972), p.70 only works for zero mean vector!
# We have to demean all variables first  
JohnsonKotzFormula <- function(mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
  lower = rep(-Inf, length = length(mean)), 
  upper = rep( Inf, length = length(mean))) {
  
  # determine which variables are truncated
  idx <- which(!is.infinite(lower) | !is.infinite(upper)) # index of truncated variables
  n <- length(mean)
  k <- length(idx) # number of truncated variables
  if (k >= n) stop(sprintf("Number of truncated variables (%s) must be lower than total number of variables (%s).", k, n))
  if (k == 0) {
   return(list(tmean=mean, tvar=sigma)) # no truncation
  }
  
  # transform to zero mean first
  lower <- lower - mean
  upper <- upper - mean
  
  # partitionining of sigma 
  # sigma = [ V11  V12 ]
  #         [ V21  V22 ]
  V11 <- sigma[idx,idx]
  V12 <- sigma[idx,-idx]
  V21 <- sigma[-idx,idx]
  V22 <- sigma[-idx,-idx]
  
  # determine truncated mean xi and truncated variance U11
  r <- mtmvnorm(mean=rep(0, k), sigma=V11, lower=lower[idx], upper=upper[idx])  
  xi <- r$tmean
  U11 <- r$tvar
  
  invV11 <- solve(V11)  # V11^(-1)
  
  # See Johnson/Kotz (1972), p.70 formula 
  tmean <- numeric(n)
  tmean[idx] <- xi
  tmean[-idx] <- xi %*% invV11 %*% V12
  tvar <- matrix(NA, n, n)
  tvar[idx, idx]  <- U11
  tvar[idx, -idx] <- U11 %*% invV11 %*% V12
  tvar[-idx, idx] <- V21 %*% invV11 %*% U11
  tvar[-idx, -idx] <- V22 - V21 %*% (invV11 - invV11 %*% U11 %*% invV11) %*% V12
  
  tmean <- tmean + mean
  
  return(list(tmean=tmean, tvar=tvar))
}

# Mean and Covariance of the truncated multivariate distribution (double truncation, general sigma, general mean)
#
# @param mean mean vector (k x 1)
# @param sigma covariance matrix (k x k)
# @param lower lower truncation point (k x 1)
# @param upper upper truncation point (k x 1)
# @param doComputeVariance flag whether to compute variance (for performance reasons)
mtmvnorm <- function(mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
  lower = rep(-Inf, length = length(mean)), 
  upper = rep( Inf, length = length(mean)),
  doComputeVariance=TRUE,
  pmvnorm.algorithm=GenzBretz())
{
  N <- length(mean)
  
  # Check input parameters
  cargs <- checkTmvArgs(mean, sigma, lower, upper)
  mean  <- cargs$mean
  sigma <- cargs$sigma
  lower <- cargs$lower
  upper <- cargs$upper
  
  # check number of truncated variables; if only a subset of variables is truncated
  # we can use the Johnson/Kotz formula together with mtmvnorm()
     
  # determine which variables are truncated
  idx <- which(!is.infinite(lower) | !is.infinite(upper)) # index of truncated variables
  k <- length(idx) # number of truncated variables
  if (k < N) {
    return(JohnsonKotzFormula(mean=mean, sigma=sigma, lower=lower, upper=upper))
  }

  # Truncated Mean
  TMEAN <- numeric(N)
  # Truncated Covariance matrix
  TVAR  <- matrix(NA, N, N)
		
  # Verschiebe die Integrationsgrenzen um -mean, damit der Mittelwert 0 wird 
  a     <- lower - mean
  b     <- upper - mean
  lower <- lower - mean
  upper <- upper - mean 
    	
  # eindimensionale Randdichte
  F_a <- numeric(N)
  F_b <- numeric(N)
  
  zero_mean <- rep(0,N)
  
  # pre-calculate one-dimensial marginals F_a[q] once
  for (q in 1:N) {
    tmp <- dtmvnorm.marginal(xn=c(a[q],b[q]), n = q, mean=zero_mean, sigma=sigma, lower=lower, upper=upper)
  	F_a[q] <- tmp[1]
 		F_b[q] <- tmp[2]
 	}
    	
  # 1. Bestimme E[X_i] = mean + Sigma %*% (F_a - F_b)
  TMEAN <- as.vector(sigma %*% (F_a - F_b))
  
  if (doComputeVariance) {
   # TODO: 
   # calculating the bivariate densities is not necessary 
   # in case of conditional independence.
   # calculate bivariate density only on first use and then cache it
   # so we can avoid this memory overhead.
   
   F2 <- matrix(0,  N, N)
   for (q in 1:N) {
     for (s in 1:N) {
       if (q != s) {
         d <- dtmvnorm.marginal2(
           xq=c(a[q], b[q], a[q], b[q]),
           xr=c(a[s], a[s], b[s], b[s]), q=q, r=s,  
           mean=zero_mean, sigma=sigma, lower=lower, upper=upper, pmvnorm.algorithm=pmvnorm.algorithm)
         F2[q,s] <- (d[1] - d[2]) - (d[3] - d[4])  
       }
     }
   }
      
   # 2. Bestimme E[X_i, X_j]
   
   # Check if a[q] = -Inf or b[q]=+Inf, then F_a[q]=F_b[q]=0, but a[q] * F_a[q] = NaN  and b[q] * F_b[q] = NaN  
   F_a_q <- ifelse(is.infinite(a), 0, a * F_a)    # n-dimensional vector q=1..N
   F_b_q <- ifelse(is.infinite(b), 0, b * F_b)    # n-dimensional vector q=1..N
  
   for (i in 1:N) {
  	for (j in 1:N) {
  	  sum <- 0
    	for (q in 1:N) {
    		sum <- sum + sigma[i,q] * sigma[j,q] * (sigma[q,q])^(-1) * (F_a_q[q] - F_b_q[q])
  			if (j != q) {
  			  sum2 <- 0
    			for (s in 1:N) {
    			    # this term tt will be zero if the partial correlation coefficient \rho_{js.q} is zero!
    				  # even for s == q will the term be zero, so we do not need s!=q condition here
    				  tt <- (sigma[j,s] - sigma[q,s] * sigma[j,q] * (sigma[q,q])^(-1))
  				    sum2 <- sum2 + tt * F2[q,s]
    			}
    			sum2 <- sigma[i, q] * sum2
    			sum <- sum + sum2
    		}	
    		} # end for q
  		TVAR[i, j] <- sigma[i, j] + sum
  		#general mean case: TVAR[i, j] = mean[j] * TMEAN[i] + mean[i] * TMEAN[j] - mean[i] * mean[j] + sigma[i, j] + sum
  	}
   }
    	
   # 3. Bestimme Varianz Cov(X_i, X_j) = E[X_i, X_j] - E[X_i]*E[X_j] für (0, sigma)-case
   TVAR <- TVAR - TMEAN %*% t(TMEAN)
  } else {
   TVAR = NA
  }
    	
  # 4. Rückverschiebung um +mean für (mu, sigma)-case
  TMEAN <- TMEAN + mean
  
  return(list(tmean=TMEAN, tvar=TVAR))
}

# Bestimmung von Erwartungswert und Kovarianzmatrix über numerische Integration und die eindimensionale Randdichte
# d.h. 
# E[X_i]       = \int_{a_i}^{b_i}{x_i * f(x_i) d{x_i}}
# Var[x_i]     = \int_{a_i}^{b_i}{(x_i-\mu_i)^2 * f(x_i) d{x_i}}
# Cov[x_i,x_j] = \int_{a_i}^{b_i}\int_{a_j}^{b_j}{(x_i-\mu_i)(x_j-\mu_j) * f(x_i,x_j) d{x_i}d{x_j}}
#
# Die Bestimmung von E[X_i] und Var[x_i]
# Die Bestimmung der Kovarianz Cov[x_i,x_j] benötigt die zweidimensionale Randdichte.
# 
#
# @param mean Mittelwertvektor (k x 1)
# @param sigma Kovarianzmatrix (k x k)
# @param lower, upper obere und untere Trunkierungspunkte (k x 1)
mtmvnorm.quadrature <- function(mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), lower = rep(-Inf, length = length(mean)), upper = rep( Inf, length = length(mean)))
{
  k       = length(mean)
  
  # Bestimmung des Erwartungswerts/Varianz über numerische Integration
  expectation <- function(x, n=1)
  {
    x * dtmvnorm.marginal(x, n=n, mean=mean, sigma=sigma, lower=lower, upper=upper)
  }

  variance <- function(x, n=1)
  {
    (x - m.integration[n])^2 * dtmvnorm.marginal(x, n=n, mean=mean, sigma=sigma, lower=lower, upper=upper)
  }

  # Determine expectation from one-dimensional marginal distribution using integration
  # i=1..k
  m.integration<-numeric(k)
  for (i in 1:k)
  {
    m.integration[i] <- integrate(expectation, lower[i], upper[i], n=i)$value 
  }
  
  # Determine variances from one-dimensional marginal distribution using integration
  # i=1..k
  v.integration<-numeric(k)
  for (i in 1:k)
  {
    v.integration[i] <- integrate(variance, lower[i], upper[i], n=i)$value 
  }
  
  return(list(m=m.integration, v=v.integration))
}
