"postmean.laplace" <-
function(x, w, a = 0.5)
{
#
# find the posterior mean for the double exponential prior for 
#   given x, w and a, assuming the error variance
#   is 1.
#
#  only allow a < 20. 
	a <- min(a, 20)	#
#  First find the odds of zero and the shrinkage factor
#
	wpost <- 1 - (1 - w)/(1 + w * beta.laplace(x, a))	
	#  now find the posterior mean conditional on being non-zero
#
	sx <- sign(x)
	x <- abs(x)
	cp1 <- pnorm(x - a)
	dp1 <- dnorm(x - a)
	cp2 <- pnorm( - x - a)
	dp2 <- dnorm(x + a)
	ef <- exp(pmin(2 * a * x, 100))
	postmeancond <- ((x - a) * cp1 + dp1 + ef * ((x + a) * cp2 - dp2))/(cp1 +
		ef * cp2)	#
#  calculate posterior mean and return
#
	mutilde <- sx * wpost * postmeancond
	return(mutilde)
}
