##' Average by session
##' 
##' Calculates the average estimate per session or block of trials
##' 
##' @param estimate	Series of estimates in event time
##' @param sessionBoundaries	Vector of the starting indices for each session
##' (which means to include the end, the last value should be length(estimate) + 1)
##' @return Vector of average estimate for each session
##' @export
##' @author Chloe Bracis
##' @examples 
##' # Create vector of values (i.e. estimates, respones, etc.)
##' values = runif(100)
##' # Specify sessions, here a group of 10 trials
##' sessionBoundaries = seq(1, 101, 10)
##' valuesBySession = averageBySession(values, sessionBoundaries)

averageBySession = function( estimate, sessionBoundaries )
{
	# can pass in NA to skip averaging by session, but still called in model fitting
	if ( all( is.na( sessionBoundaries ) ) )
	{
		avg = estimate
	} else
	{
		avg = sapply( 1:( length( sessionBoundaries ) - 1 ),
					  function( x ) if ( !is.na( sessionBoundaries[x] ) 
					  				   && !is.na( sessionBoundaries[x + 1] ) )
					  	mean( estimate[sessionBoundaries[x]:( sessionBoundaries[x + 1] - 1 )] ) )
	}
	return( unlist( avg ) )
}
