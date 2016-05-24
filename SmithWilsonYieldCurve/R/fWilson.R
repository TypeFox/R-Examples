#' Wilson function
#' 
#' Acts as a kernel for regression
#' 
#' @param t a time
#' @param u another time
#' @param ufr the ultimate forward rate
#' @param alpha the spead of reversion to the ultimate forward rate
#' 
fWilson <- function(t, u, ufr, alpha) {

	min_tu <- outer( t, u, pmin )
	max_tu <- outer( t, u, pmax)
	
	fw <- function( t, u ) exp( -ufr*(t+u) ) * ( alpha*min_tu - 0.5 * ( exp( alpha*min_tu ) - exp( -alpha*min_tu) ) / exp( alpha*max_tu ) )
	
	W <- t( outer( t, u, fw ) )
	
	return( W )
	
}