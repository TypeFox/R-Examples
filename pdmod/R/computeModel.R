##' Calculates proximal/distal model
##' 
##' Calulates a realization of a proximal/distal model for a specified sequence of trials and paramter values. Use the \code{verbose}
##' parameter to include underlying model components (distal and proximal estimates, weights, uncertainties and signal-reward 
##' association) in addition to the mean estimate.
##' 
##' @param x			Object of class \code{\link{TimedVector}} specifying trials including 
##' whether signal was rewarded/unrewarded and times
##' @param mFast		Learning rate of proximal memory estimates
##' @param mSlow		Learning rate of distal memory estimates
##' @param n			Learning rate of uncertainty estimates
##' @param h			Decay rate of distal memory uncertainty estimator as time passes between trials
##' @param g			Association learning speed parameter
##' @param tau		Temporal scaling coefficient to translate time differences in \code{x} to fractional days.
##' Defaults to \code{1/TV_DAY} assuming that the times in \code{x} are expressed in minutes.
##' @param threshold	Difference in real time that must pass before deflation kicks in (used for testing)
##' @param verbose	true to include supporting estimates, weights, etc.
##' @return Series of estimates
##' @seealso \code{\link{calculateResponse}}, \code{\link{averageBySession}}
##' @export
##' @useDynLib pdmod compute_model
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
 
computeModel = function( x, mFast, mSlow, n, g = 0, 
 						 h, tau = 1/TV_DAY, threshold = 0, verbose = TRUE  )
{
	if ( !isTimedVector( x ) )
		stop( "x must be of class TimedVector" )
	if ( mFast < 0 || mFast > 1 )
		stop( "mFast must be between 0 and 1" )
	if ( mSlow < 0 || mSlow > 1 )
		stop( "mSlow must be between 0 and 1" )
	if ( n < 0 || n > 1 )
		stop( "n must be between 0 and 1" )
	if ( h < 0 || h > 1 || !is.numeric( h ) )
		stop( "h must be between 0 and 1" )
	if ( g < 0 || !is.numeric( g ) )
		stop( "g must be greater than or equal to 0" )
	
	len = length( x )
	
	interference = rep( 0, len ) 
	a = 1
	b = 0
	t = time( x )
	deltat = vector( "numeric", length( t ) )
	deltat[2:length( t )] = t[2:length( t )] - t[1:( length( t ) - 1 )]
	deltat[1] = t[1]
	
	result = .C( "compute_model",
				 as.double( x ),
				 as.double( deltat ),
				 as.integer( interference ),
				 as.integer( len ),
				 as.double( mFast ),
				 as.double( mSlow ),
				 as.double( n ),
				 as.double( g ),
				 as.double( h ),
				 as.double( tau ),
				 as.double( threshold ),
				 as.double( a ),
				 as.double( b ),
				 est = double( 2 * len ),
				 estMse = double( 2 * len ),
				 weights = double( 2 * len ),
				 forecastMse = double( len ),
				 y = double( len ),
				 forecast = double( len ),
				 NAOK = TRUE )
	
	forecast = TimedVector( result$forecast, time( x ) )
	
	if ( verbose )
	{
		attr( forecast, "class") = c( "TimedVector", "pdmod" )
		attr( forecast, "estimates" ) = matrix( result$est, ncol = 2 )
		attr( forecast, "errorEstimates" ) = matrix( result$estMse, ncol = 2 )
		attr( forecast, "weights" ) = matrix( result$weights, ncol = 2 )
		attr( forecast, "y" ) = result$y
		attr( forecast, "forecastErrorEstimates" ) = result$forecastMse
	}
	
	return( forecast )
}
