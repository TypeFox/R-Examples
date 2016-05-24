##' Create a TimedVector
##' 
##' The class TimedVector contains a vector of values in event time, 
##' as well as when in real time those events took place.
##' 
##' @param x Series of values in event time
##' @param t (optional) Cooresponding real time of events in minutes.  Default is an event every minute.
##' @return TimedVector
##' @export
##' @seealso \code{\link{Constants}}, \code{\link{isTimedVector}}, \code{\link{verifyTimedVector}}
##' @author Chloe Bracis
##' @examples
##' # One session of 20 rewarded trials every minute
##' TimedVector(rep(1, 20), 1:20)
##' 
##' # Three sessions of rewarded trials, then one session of non-rewarded trials, 
##' # with trials every 2 min and sessions every day
##' trialTime = as.vector(sapply(0:3, function(x) seq(2, 20, 2) + x * TV_DAY))
##' TimedVector(c(rep(1, 30), rep(0, 10)), trialTime)
##' 
##' # The above schedule of sessions, but 50% probability of reward
##' TimedVector(sample(0:1, 40, replace = TRUE), trialTime)

TimedVector = function( x, t )
{
  if ( is.null( x ) && is.null( t ) )
    return( NULL )
  if ( missing( t ) )
    t = 1:length( x )
  stopifnot( length( x ) == length( unique( t ) ) )
  stopifnot( all ( t == sort( t ) ) )
  
  tv = x
  attr( tv, "class" ) = "TimedVector"
  attr( tv, "realTime" ) = t
  return( tv )
}
