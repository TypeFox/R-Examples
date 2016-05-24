formula.systemfit <- function( x, ... ) {
   result <- list()
   eqnLabels <- NULL
   for( i in 1:length( x$eq ) ){
      result <- c( result, formula( x$eq[[ i ]] ) )
      eqnLabels <- c( eqnLabels, x$eq[[ i ]]$eqnLabel )
   }
   names( result ) <- eqnLabels
   return( result )
}

formula.systemfit.equation <- function( x, ... ) {
   result <- formula( x$terms )
   return( result )
}
