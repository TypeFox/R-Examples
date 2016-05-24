summary.aidsElas <- function( object, ... ) {

   result <- object
   allCoef <- c( object$exp, t( object$hicks ), t( object$marshall ) )

   if( is.null( object$allVcov ) ) {
      result$table <- matrix( NA, length( allCoef ), 4 )
      result$table[ , 1 ] <- allCoef
      colnames( result$table ) <-
         c( "Estimate", "Std. Error", "t value", "Pr(>|t|)" )
   } else {
      result$table <- coefTable( allCoef,
         diag( object$allVcov )^0.5, result$df )
      rownames( result$table ) <- rownames( object$allVcov )
   }

   class( result ) <- "summary.aidsElas"
   return( result )
}