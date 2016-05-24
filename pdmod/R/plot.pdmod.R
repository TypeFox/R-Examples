##' Plot model
##' 
##' Plots a \code{pdmod} class (what's returned from \code{\link{computeModel}} with \code{verbose = TRUE}). 
##' The plots show the proximal and distal estimates, their corresponding uncertainties and weights, 
##' as well as the overall mean estimate.
##' 
##' @rdname plot.pdmod
##' @method plot pdmod
##' @S3method plot pdmod
##' @export
##' @param x		Object of class \code{pdmod}
##' @param actual	Actual rewards received
##' @param n		(optional) Only plot the last n values
##' @param ...		Other arguments to \code{\link{plot}}
##' @author Chloe Bracis
##' @examples 
##' # Create 5 sessions of 20 rewarded trials, 
##' # then 2 sessions of 20 unrewarded trials
##' trialTime = as.vector(sapply(0:6, function(x) 1:20 + x * TV_DAY))
##' trials =  TimedVector(c(rep(1, 5*20), rep(0, 2*20)), trialTime)
##' 
##' estimates = computeModel(trials, mFast = 0.7, mSlow = 0.1, n = 0.05, 
##' 						 g = 500, h = 0.2, verbose = TRUE)
##' plot(estimates, trials)

plot.pdmod = function( x, actual, n, ... )
{
	stopifnot(inherits( x, "pdmod" ))
	est = attr( x, "estimates" )
	mse = attr( x, "errorEstimates" )
	weights = attr( x, "weights" )
	total = length( est[,1] )
	if ( missing( n ) )
		n = total
	xlim = c( total - n + 1, total )	
	xs = ( total - n + 1 ):total
	
	
	oldPar = par( mfrow = c( 2, 2 ), mar = c( 3, 3, 2, 1 ) + 0.1, mgp = c( 1.5, 0.5, 0 ) )
	plot(0, 0, type = "n", ylim = c( 0, max( est, na.rm = TRUE ) ), xlim = xlim, 
		 ylab = "Estimate", xlab = "", main = "Estimates", ... )
	lines( xs, tail( est[,1], n ), lty = 1, col = 3 )
	lines( xs, tail( est[,2], n ), lty = 1, col = 2 )
	legend( "topleft",  bty = "n",
			legend = c( "Proximal", "Distal" ),
			lty = 1, col = c( 3, 2 ) )
	
	plot(0, 0, type = "n", ylim = c( 0, 1 ), xlim = xlim, 
		 ylab = "Weight", xlab = "", main = "Weights", ... )
	lines( xs, tail( weights[, 1], n ), lty = 1, col = 3 )
	lines( xs, tail( weights[, 2], n ), lty = 1, col = 2 )
	
	plot(0, 0, type = "n", ylim = c( 0, max( mse, na.rm = TRUE ) ), xlim = xlim, 
		 ylab = "Uncertainty", xlab = "Trial", main = "Uncertainties", ... )
	lines( xs, tail( mse[,1], n ), lty = 1, col = 3 )
	lines( xs, tail( mse[,2], n ), lty = 1, col = 2 )
	
	plot(0, 0, type = "n", ylim = c( 0, max( max( x, na.rm = TRUE ), max( actual, na.rm = TRUE ) ) ), xlim = xlim, 
		 ylab = "Mean estimate", xlab = "Trial", main = "Mean estimate", ... )
	lines( xs, tail( x, n ), lty = 1 )
	points( xs, tail( actual, n ) )
	legend( "topleft", bty = "n",
			legend = c( "Mean estimate", "Observed" ),
			lty = c( 1, NA ), pch = c( NA, 1 ) )
	par( oldPar )
}
