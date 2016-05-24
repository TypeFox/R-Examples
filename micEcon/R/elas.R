## Return elasticities

elas <- function( object, ... )
    UseMethod("elas")

elasticities <- function( object, ... )
    UseMethod("elas")

elas.default <- function( object, ... ) {
   return( object$elas )
}