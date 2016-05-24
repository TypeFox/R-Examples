fitted.probit <- function( object, ... ) {
   result <- drop( pnorm( linearPredictors( object ) ) )
   names( result ) <- row.names( model.matrix( object ) )
   return( result )
}
