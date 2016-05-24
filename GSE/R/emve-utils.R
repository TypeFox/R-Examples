.emve.Rcpp <- function(x, x_nonmiss, pu, n, p, theta0, G, d, x.miss.group.match, miss.group.unique, miss.group.counts,
		miss.group.obs.col, miss.group.mis.col, miss.group.p, miss.group.n, nresample, nsubsize, minrcond, cc, ck, EM_maxits){

	res <- tryCatch( .Call("emve_Rcpp", x, x_nonmiss, pu, n, p, theta0, G, d, x.miss.group.match, miss.group.unique, miss.group.counts,
		miss.group.obs.col, miss.group.mis.col, miss.group.p, miss.group.n, nresample, nsubsize, minrcond, cc, ck, EM_maxits),
		"std::range_error" = function(e){
		conditionMessage( e ) } )

	if( is.character( res ) )
		stop( paste(res, "\nEstimated covariance matrix is not positive definite.\n") )

	ss0 <- res[1,1,]
	mu0 <- t(res[2,,])
	S0 <- res[-c(1,2),,]

	list(ss0=ss0, mu0=mu0, S0=S0)
}

###################################################################
## Compute the correction constants k and 50% chi-sq quantiles for 
## the computation of the scale of MVE for incomplete data (new version)
.scale.mve.init <- function(pp){
	cc <- qchisq(1/2,pp)
	k1 <- 1/pp
	k2 <- 0.5^(pp/2)
	k3 <- 1/(gamma(pp/2))
	k4 <- cc^(1+pp/2)
	k5 <- exp(-cc/2)
	kk <- k1 * k2 * k3 * k4 * k5
	ck <- cc*kk
	return( list(cc=cc, ck=ck) )
}


###################################################################
## Partial mahalanobis distance
## Rcpp version
.partial.mahalanobis.Rcpp <- function(x_mu_diff, Sigma, miss_group_unique, miss_group_counts){
	res <- tryCatch( .Call("fast_partial_mahalanobis", x_mu_diff, Sigma, miss_group_unique, miss_group_counts),
		"std::range_error" = function(e){
		conditionMessage( e ) } )
	if( is.character( res ) ) stop(res)
	return( c(res) )
}






	
	
	