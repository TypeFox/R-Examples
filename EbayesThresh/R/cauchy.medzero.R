"cauchy.medzero" <-
function(x, z, w)
{
# the objective function that has to be zeroed, component by component, 
# to find the 
#  posterior median when the quasi-Cauchy prior is used. 
#   x is the parameter vector, z is the data vector, w
# is the weight
#   x and z may be scalars
#
	hh <- z - x
	dnhh <- dnorm(hh)
	yleft <- pnorm(hh) - z * dnhh + ((z * x - 1) * dnhh * pnorm( - x))/
		dnorm(x)
	yright2 <- 1 + exp( - z^2/2) * (z^2 * (1/w - 1) - 1)
	return(yright2/2 - yleft)
}
