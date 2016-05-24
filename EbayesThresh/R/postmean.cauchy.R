"postmean.cauchy" <-
function(x, w)
{
#
#  Find the posterior mean for the quasi-Cauchy prior with mixing weight w
#   given data x, which may be a scalar or a vector.
#
	muhat <- x
	ind <- (x == 0)
	x <- x[!ind]
	ex <- exp( - x^2/2)
	z <- w * (x - (2 * (1 - ex))/x)
	z <- z/(w * (1 - ex) + (1 - w) * ex * x^2)
	muhat[!ind] <- z
	return(muhat)
}
