##' Fit model parameters
##' 
##' Estimates parameters for proximal/distal model using multi-criteria estimation (\code{\link[=nsga2]{mco}})
##' 
##' @param dataX			Object of class \code{\link{TimedVector}} specifying trials including 
##' whether signal was rewarded/unrewarded and times
##' @param dataResponse		Corresponding observations of subject's response to signal
##' @param responseFunction	The function to use to transform the mean estimate into a response
##' @param sessionBoundaries	(optional) Vector defining how to group the trials into sessions
##' where the items are the starting indicies for each session (so the last value can be the index 
##' after the last trial) and \code{NA}s are used for gaps between sessions
##' @param fitG				\code{TRUE} (default) to estimate g, or \code{FALSE} to fix g at 0
##' @return Model fit
##' @seealso \code{\link{computeModel}}
##' @export
##' @importFrom mco nsga2
##' @author Chloe Bracis
fitModel = function( dataX, dataResponse, 
						responseFunction = calculateResponse, 
						sessionBoundaries = NA, fitG = TRUE )
{
	nDatasets = length( dataX )
	nParameters = 6
	if ( fitG ) { nParameters = nParameters + 1 }
	upperBounds = c( 1, 1, 1, 1, 10, 100 )
	if ( fitG ) { upperBounds = c( upperBounds, 10000 ) }
	lowerBounds = rep( 0, nParameters )
	
	fit = nsga2( modelObjectiveFunction, idim = nParameters, odim = nDatasets, 
				 nDatasets, dataX, dataResponse, responseFunction, sessionBoundaries, fitG,
				 lower.bounds = lowerBounds, 
				 upper.bounds = upperBounds,
				 popsize = 200, generations = 200 )
	return( fit )
}
