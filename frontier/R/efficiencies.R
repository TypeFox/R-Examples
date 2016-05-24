## Return efficiency estimates

efficiencies <- function( object, ... )
    UseMethod( "efficiencies" )

# default method
efficiencies.default <- function( object, ... ) {
   return( object$effic )
}
