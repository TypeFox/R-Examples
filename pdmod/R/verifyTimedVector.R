##' Verify TimedVector
##' 
##' Verifies object really is a \code{\link{TimedVector}} (stronger checks than \code{\link{isTimedVector}}).
##' 
##' @param tv Potential \code{TimedVector} object
##' @return \code{TRUE}, if the object is a \code{TimedVector}
##'			\code{FALSE}, otherwise
##' @seealso \code{\link{isTimedVector}}, \code{\link{TimedVector}}
##' @export
##' @author Chloe Bracis
verifyTimedVector = function( tv )
{
	if ( !isTimedVector( tv ) ) 
		return( FALSE )
	if ( is.null( attr( tv, "realTime" ) ) )
		return( FALSE )
	
	t = attr( tv, "realTime" )
	if ( length( tv ) != length( unique( t ) ) )
		return( FALSE )	
	if ( !all( t == sort( t ) ) )
		return( FALSE )
	
	return( TRUE )
}
