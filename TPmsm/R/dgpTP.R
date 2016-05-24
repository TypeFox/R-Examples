dgpTP <- function(n, corr, dist, dist.par, model.cens, cens.par, state2.prob) {
	if ( missing(n) ) stop("Argument 'n' is missing, with no default")
	if ( missing(corr) ) stop("Argument 'corr' is missing, with no default")
	if ( missing(dist) ) stop("Argument 'dist' is missing, with no default")
	if ( missing(dist.par) ) stop("Argument 'dist.par' is missing, with no default")
	if ( missing(model.cens) ) stop("Argument 'model.cens' is missing, with no default")
	if ( missing(cens.par) ) stop("Argument 'cens.par' is missing, with no default")
	if ( missing(state2.prob) ) stop("Argument 'state2.prob' is missing, with no default")
	return( .Call(Rf_dgpTP, as.integer(n), as.double(corr), dist, as.double(dist.par), model.cens, as.double(cens.par), as.double(state2.prob), PACKAGE="TPmsm") )
}
