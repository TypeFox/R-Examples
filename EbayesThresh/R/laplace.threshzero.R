"laplace.threshzero" <-
function(x, w, a = 0.5)
{
#
# the function that has to be zeroed to find the threshold with the Laplace
#    prior.  
#  only allow a < 20. 
	a <- min(a, 20)	#
	z <- pnorm(x - a) - (dnorm(x - a) * (1/w + beta.laplace(x, a)))/a
	return(z)
}
