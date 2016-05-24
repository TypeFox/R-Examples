source("R/rtmvnorm.R")
# Dichtefunktion der Multivariaten Trunkierten Normalverteilung mit Trunkierungsvektor lower and upper
#
# vgl. Horrace (2005) "Some Results on the Multivariate Truncated Normal Distribution"
#
# @param x Argumentenvektor der Dichte der Länge n oder Matrix (T x n) mit T Beobachtungen
# @param mean  Mittelwertvektor der Länge n
# @param sigma Kovarianzmatrix (n x n)
# @param lower unterer Trunkierungsvektor (n x 1) mit lower <= x <= upper
# @param upper oberer Trunkierungsvektor (n x 1) mit lower <= x <= upper
# @param margin if NULL then joint density, if MARGIN=1 then first marginal density, if MARGIN=c(1,2) 
#               then bivariate marginal density for x_1 and x_2
dtmvnorm <- function(x, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
  lower = rep( -Inf, length = length(mean)), upper = rep( Inf, length = length(mean)), 
		 log = FALSE, margin=NULL)
{
  # check of standard tmvnorm arguments
  cargs <- checkTmvArgs(mean=mean, sigma=sigma, lower=lower, upper=upper)
  mean  <- cargs$mean
  sigma <- cargs$sigma
  lower <- cargs$lower
  upper <- cargs$upper
  
  # Check of optional argument "margin"
  if (!is.null(margin)) {
	# Aufpassen! dtmvnorm() nimmt als Argumente auch eine (T x n)-Matrix,
	# dtmvnorm.marginal() nimmt nur einen Vektor
	# dtmvnorm.marginal2() nimmt 2 Vektoren der gleichen Länge
	# Aufpassen mit Checks auf die Länge von x
	# Aufpassen mit dem log=TRUE Argument!
    if (!length(margin) %in% c(1, 2))
	  stop("Length of margin must be either 1 (one-dimensional marginal density) or 2 (bivariate marginal density).")
	if (any(margin <= 0) || any(margin > length(mean))) {
	  stop("All elements in margin must be in 1..length(mean).")	
	}
	# one-dimensional marginal density f_{n}(x_n)
	if (length(margin) == 1) {
	  return(dtmvnorm.marginal(xn=x, n=margin, mean = mean, sigma = sigma, lower = lower, upper = upper, log = log))		
	}
	# for bivariate marginal density f_{q,r}(x_q, x_r) we need q <> r and "x" as (n x 2) matrix
	if (length(margin) == 2) {
	  if(margin[1] == margin[2])	
	    stop("Two different margins needed for bivariate marginal density.")
	  if (is.vector(x)) {
	    x <- matrix(x, ncol = length(x))
    }  
	  if(!is.matrix(x) || ncol(x) != 2)
		stop("For bivariate marginal density x must be either a (n x 2) matrix or a vector of length 2.")  
	  # bivariate marginal density f_{q,r}(x_q, x_r)	
	  return(dtmvnorm.marginal2(xq=x[,1], xr=x[,2], q=margin[1], r=margin[2], mean = mean, sigma = sigma, lower = lower, upper = upper, log = log))	  
    }	
  }
  
  # Check of additional inputs like x
  if (is.vector(x)) {
	  x <- matrix(x, ncol = length(x))
  }
    
  # Anzahl der Beobachtungen
  T <- nrow(x)
  
  # check for each row if in support region
  insidesupportregion <- logical(T)
  for (i in 1:T)
  {
    insidesupportregion[i] = all(x[i,] >= lower & x[i,] <= upper & !any(is.infinite(x)))
  }
  
  if(log) {
    # density value for points inside the support region
    dvin <- dmvnorm(x, mean=mean, sigma=sigma, log=TRUE) - log(pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma)) 
    # density value for points outside the support region
    dvout <- -Inf
  } else {
    dvin <- dmvnorm(x, mean=mean, sigma=sigma, log=FALSE) / pmvnorm(lower=lower, upper=upper, mean=mean, sigma=sigma)
    dvout <- 0
  }
  
  
  f <- ifelse(insidesupportregion, dvin, dvout)
  return(f)
}

#dtmvnorm(x=c(0,0))
#dtmvnorm(x=c(0,0), sigma=diag(2))
#dtmvnorm(x=c(0,0), mean=c(0,0), sigma=diag(2))
#dmvnorm(x=c(0,0), mean=c(0,0), sigma=diag(2))
#dtmvnorm(x=matrix(c(0,0,1,1),2,2, byrow=TRUE), mean=c(0,0), sigma=diag(2))
#dtmvnorm(x=matrix(c(0,0,1,1),2,2, byrow=TRUE), mean=c(0,0), sigma=diag(2), lower=c(-1,-1), upper=c(0.5, 0.5))
#dtmvnorm(x=matrix(c(0,0,1,1),2,2, byrow=TRUE), mean=c(0,0), sigma=diag(2), lower=c(-1,-1), upper=c(0.5, 0.5), log=TRUE)
#dtmvnorm(as.matrix(seq(-1,2, by=0.1), ncol=1), mean=c(0.5), sigma=as.matrix(1.2^2), lower=0)
