##' @method c TimedVector
##' @S3method c TimedVector
c.TimedVector = function( ... )
{
	args = list(...)
	if( !all( sapply( args, isTimedVector ) ) )
	{
		return( unlist( sapply( args, as.vector ) ) )
	}
	totalLength = sum( sapply( args, length ) )
	x = vector( "numeric", length = totalLength )
	t = vector( "numeric", length = totalLength )
	
	index = 1
	for ( i in 1:length( args ) )
	{
		len = length( args[[i]] )
		x[index:(index + len - 1)] = args[[i]]
		t[index:(index + len - 1)] = attr( args[[i]], "realTime" )
		index = index + len
	}
	
	tv = TimedVector( x, t )
	return( tv )
}
