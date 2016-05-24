##' @method print TimedVector
##' @S3method print TimedVector
print.TimedVector = function( x, ... )
{
	print( cbind( time = attr( x, "realTime" ), value = x ), ... )
}
