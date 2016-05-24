terms.systemfit <- function( x, ... ) {
   result <- list()
   eqnLabels <- NULL
   for( i in 1:length( x$eq ) ){
      result <- c( result, terms( x$eq[[ i ]] ) )
      eqnLabels <- c( eqnLabels, x$eq[[ i ]]$eqnLabel )
   }
   names( result ) <- eqnLabels
   return( result )
}

terms.systemfit.equation <- function( x, ... ) {
   result <- x$terms
   return( result )
}
