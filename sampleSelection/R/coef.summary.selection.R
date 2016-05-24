coef.summary.selection <- function( object, part="full", ... ) {

   if( !( part %in% c( "full", "outcome" ) ) ) {
      stop( "argument 'part' must be either 'full' or 'outcome'" )
   }

   result <- object$estimate

   if( part == "outcome" ) {
      result <- result[ object$param$index$outcome, ]
   }

   return( result )
}
