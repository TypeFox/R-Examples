dgBIV <- function(n, dist, corr, dist.par) {
	return( .Call("dgBIV", as.integer(n), dist, as.double(corr), as.double(dist.par), PACKAGE="genSurv") )
}
