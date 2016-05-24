# Density function for the truncated multivariate t-distribution
# 
# Author: stefan
###############################################################################

# Density function for the truncated multivariate t-distribution
# @param x
# @param mean
# @param sigma
# @param df degrees of freedom parameter
# @param log
dtmvt <- function(x, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), df = 1, lower= rep( -Inf, length = length(mean)), 
		upper = rep( Inf, length = length(mean)), log = FALSE){
	
	# Check of additional inputs like x
	if (is.vector(x)) {
		x <- matrix(x, ncol = length(x))
	}
	
	# Anzahl der Beobachtungen
	T = nrow(x)
	
	# check for each row if in support region
	insidesupportregion <- logical(T)
	for (i in 1:T)
	{
		insidesupportregion[i] = all(x[i,] >= lower & x[i,] <= upper & !any(is.infinite(x)))
	}
	
	# density value for points outside the support region
	dv = if (log) { -Inf } else { 0 }
	
	# conditional density
	f <- ifelse(insidesupportregion, 
	  dmvt(x, delta=mean, sigma=sigma, df=df, log=log) / pmvt(lower=lower, upper=upper, delta=mean, sigma=sigma, df=df, type="shifted"), dv)
	return(f)
}

if (FALSE) {
# Example
x1<-seq(-2, 3, by=0.1)
x2<-seq(-2, 3, by=0.1)

mean=c(0,0)
sigma=matrix(c(1, -0.5, -0.5, 1), 2, 2)
lower=c(-1,-1)


density<-function(x)
{
	z=dtmvt(x, mean=mean, sigma=sigma, lower=lower)
	z
}

fgrid <- function(x, y, f)
{
	z <- matrix(nrow=length(x), ncol=length(y))
	for(m in 1:length(x)){
		for(n in 1:length(y)){
			z[m,n] <- f(c(x[m], y[n]))
		}
	}
	z
}

# compute multivariate-t density d for grid
d=fgrid(x1, x2, function(x) dtmvt(x, mean=mean, sigma=sigma, lower=lower))

# compute multivariate normal density d for grid
d2=fgrid(x1, x2, function(x) dtmvnorm(x, mean=mean, sigma=sigma, lower=lower))

# plot density as contourplot
contour(x1, x2, d, nlevels=5, main="Truncated Multivariate t Density", 
		xlab=expression(x[1]), ylab=expression(x[2]))

contour(x1, x2, d2, nlevels=5, add=TRUE, col="red")
abline(v=-1, lty=3, lwd=2)
abline(h=-1, lty=3, lwd=2)
}
