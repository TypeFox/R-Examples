setPackageSeedTP <- function(seed=12345) {
	if ( !is.numeric(seed) ) stop("Argument 'seed' must be numeric")
	ll <- length(seed)
	if (ll < 6) {
		tmpseed<-12345
		seed <- c( seed, rep(tmpseed,6-ll) )
	}
	if (ll > 6) {
		seed <- seed[1:6]
		warning( "Seed may not exceed the length of 6. Truncated to ", paste(seed, collapse=", ") )
	}
	return( invisible( .Call(Rf_rset_package_seed, as.double(seed), PACKAGE="TPmsm") ) )
}
