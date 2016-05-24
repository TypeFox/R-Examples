r.gaussmix <-
function(n, means, prior = rep(1/length(means), length(means)))
{
	###density function for a standardized mixture of gaussians
	xbar <- sum(means * prior)
	xvar <- 1 + sum(prior * (means - xbar)^2)
	xsd <- sqrt(xvar)
	nj <- drop(rmultinom(1,n, matrix(prior, 1, length(prior))))
	x <- NULL
	for(i in seq(prior)) {
		x <- c(x, rnorm(nj[i]) + means[i])
	}
	sample((x - xbar)/xsd)
}

