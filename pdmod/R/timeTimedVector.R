##' @method time TimedVector
##' @S3method time TimedVector
##' @importFrom stats time

time.TimedVector = function( x, ... )	
{
	
	return( attr( x, "realTime" ) )
	
}
