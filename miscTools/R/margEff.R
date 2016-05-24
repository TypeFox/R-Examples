margEff <- function( object, ... )
    UseMethod( "margEff" )

# default method
margEff.default <- function( object, ... ) {
   stop( "there is currently no default method available" )
}

