##' Is TimedVector
##' 
##' Determines if an object inherits \code{\link{TimedVector}}
##' 
##' @param tv Potential \code{TimedVector} object
##' @return \code{TRUE}, if the object inherits \code{TimedVector}
##'			\code{FALSE}, otherwise
##'	@seealso \code{\link{TimedVector}}, \code{\link{verifyTimedVector}}
##' @export
##' @author Chloe Bracis
##' @examples
##' # A TimedVector
##' tv = TimedVector(rep(1, 10), 1:10)
##' isTimedVector(tv)
##' 
##' # Not a TimedVector
##' isTimedVector(1:10)
##' isTimedVector(time(tv))

isTimedVector = function( tv )
{
	if ( is.null( tv ) ) 
		return( FALSE )
	if ( !inherits( tv, "TimedVector" ) ) 
		return( FALSE )
	
	return( TRUE )
}
