## create the summary results
summary.solveLP <- function(object,...) {
   class( object ) <- c( "summary.solveLP", class( object ) )
   return( object )
}
