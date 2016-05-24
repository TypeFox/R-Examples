# Berechnet die Quantile der eindimensionalen Randverteilung über uniroot()
#
# @param p probability
# @param interval a vector containing the end-points of the interval to be searched by uniroot.
# @param tail specifies which quantiles should be computed. lower.tail gives the quantile x for which P[X <= x] = p, upper.tail gives x with P[X > x] = p and both.tails leads to x with P[-x <= X <= x] = p.
# @param n
# @param mean
# @param sigma
# @param lower
# @param upper
# @param ... additional parameters to uniroot()
qtmvnorm.marginal <- function (p, 
		interval = c(-10, 10), 
		tail = c("lower.tail", "upper.tail", "both.tails"), 
		n=1, 
		mean=rep(0, nrow(sigma)), 
		sigma=diag(length(mean)), 
		lower=rep(-Inf, length = length(mean)), 
		upper=rep( Inf, length = length(mean)),
		...)
{
	if (length(p) != 1 || (p <= 0 || p >= 1)) 
		stop(sQuote("p"), " is not a double between zero and one")
	
	if (n > length(mean) || n < 1) 
		stop(sQuote("n"), " is not a integer between 1 and ",length(mean))
	
	pfct <- function(q) 
	{
		switch(tail, both.tails = {
					low <- lower
					low[n] <- -abs(q)
					upp <- upper
					upp[n] <- abs(q)
				}, upper.tail = {
					low <- lower
					upp <- upper
					low[n] <- q
				}, lower.tail = {
					low <- lower
					upp <- upper
					upp[n] <- q
				}, )	
		ptmvnorm(low, upp, mean, sigma, lower, upper) - p	
	}
	qroot <- uniroot(pfct, interval = interval, ...)
	qroot
}	