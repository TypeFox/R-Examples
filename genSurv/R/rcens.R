runifcens <- function(n, cens.par) {
	return( runif(n=n, min=0, max=cens.par) )
}

rexpcens <- function(n, cens.par) {
	return( rexp(n=n, rate=1/cens.par) )
}
