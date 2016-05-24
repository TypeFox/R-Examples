##' Objective function to fit model parameters
##' 
##' Function passed to optimization routine to minimize to estimate parameters.  
##' Uses mean squared error to  calculate difference between \code{dataResponse} and 
##' what \code{\link{computeModel}}) would forcast for \code{dataX} using parameters \code{pars}.
##' 
##' @param pars				Vector of parameters mFast, mSlow, n, hSlow, and r
##' @param dimension		What dimension to return error in, 1 for single
##' criteria optimization, or number of columns of data for multicriteria optimization
##' @param dataX			List of observations of process x(i) (with real time)
##' @param dataResponse		Corresponding list of observations of subject's response
##' to x(i), i.e. ~x(i)
##' @param responseFunction	The function to use to transform the forecast into a response
##' @param sessionBoundaries	(option) Vector defining how to group the trials into sessions
##' where the items are the starting indicies for each session (so the last value can be the 
##' index after the last trial) and \code{NA}s are
##' used for gaps between sessions
##' @param fitG				\code{TRUE} to estimate g, or \code{FALSE} to fix g at 0
##' @return Error between \code{dataRespones} and what would have been
##' estimated for \code{dataX} based on parameters pars
##' @seealso \code{\link{computeModel}}, \code{\link{fitModel}}
##' @export
##' @author Chloe Bracis
modelObjectiveFunction = function( pars, dimension, dataX, dataResponse, 
								   responseFunction = calculateResponse, sessionBoundaries = NA, fitG = TRUE )
{
	mFast = pars[1]
	mSlow = pars[2]
	n = pars[3]
	h = pars[4]
	tau = 1/TV_DAY # fix tau to convert minutes to days
	r = pars[5]
	rmax = pars[6]
	g = if ( fitG ) pars[7]	else 500
	
	if ( !is.list( dataX ) || !is.list( dataResponse ) )
		stop( "data and responses must be passed as lists" )
	
	nData = length( dataX )
	
	if ( dimension != 1 && dimension != nData )
		stop( "dimension should be 1 or the number of data sets" )
	badParameter = rep( Inf, dimension )
	
	if ( mFast < 0 || mFast > 1 || mSlow < 0 || mSlow > 1 )
		return( badParameter )
	if( mFast < mSlow )
		return( badParameter )
	if ( n < 0 || n > 1 )
		return( badParameter )
	if ( h < 0 || h > 1 )
		return( badParameter )
	#	if ( tau < 0 || tau > 1 )
	#		return( badParameter )
	if ( is.nan( n ) || is.nan( g ) || is.nan( h ) || is.nan( r ) || is.nan( rmax )  )
		return( badParameter )
	
	est = lapply( 1:nData, function ( x ) 
		computeModel( x = dataX[[x]], 
					  mFast = mFast, 
					  mSlow = mSlow,
					  n = n, 
					  g = g,
					  h = h,
					  tau = tau,
					  verbose = FALSE ) )
	
	
	# average across trials if data is reported in blocks/sessions, first calculating response from estimate
	response = lapply( 1:nData, function( x ) responseFunction( r, rmax, est[[x]] ) )
	
	if ( is.list( sessionBoundaries ) )
		response = lapply( 1:nData, function( x ) averageBySession( response[[x]], sessionBoundaries[[x]] ) )
	else
		response = lapply( 1:nData, function( x ) averageBySession( response[[x]], sessionBoundaries ) )
	
	mse = sapply( 1:nData, function ( x ) mean( ( response[[x]] - dataResponse[[x]] )^2 ) )
	if ( dimension == 1 )
		mse = sum( mse )
	
	return( mse )
}
