d.gaussmix <-
function(x, means, prior = rep(1/length(means), length(means)))
{
###density function for a standardized mixture of gaussians
	xbar <- sum(means * prior)
	xvar <- 1. + sum(prior * (means - xbar)^2.)
	xsd <- sqrt(xvar)
	dx <- x * 0.
	for(i in seq(prior)) {
		dx <- dx + prior[i] * dnorm(x * xsd - means[i] + xbar) * xsd
	}
	dx
}

