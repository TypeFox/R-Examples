"cauchy.threshzero" <-
function(z, w)
{
# the objective function that has to be zeroed
# to find the Cauchy
#  threshold.  z is the putative threshold vector, w
# is the weight
#   w can be a vector
#
	y <- pnorm(z) - z * dnorm(z) - 1/2 - (z^2 * exp( - z^2/2) * (1/w - 1))/
		2
	return(y)
}
