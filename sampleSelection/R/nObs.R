nObs.probit <- function(x, ...) {
   return( x$param$nObs )
}

nobs.probit <- function( object, ... ) {
   return( nObs.probit( object, ... ) )
}

nObs.selection <- function(x, ...) {
   return( x$param$nObs )
}

nobs.selection <- function( object, ... ) {
   return( nObs.selection( object, ... ) )
}
