logLik.mvProbit <- function( object, coef = NULL, data = NULL, 
      algorithm = NULL, nGHK = NULL, random.seed = NULL, ... ) {

   if( is.null( coef ) && is.null( data ) &&
         is.null( algorithm ) && is.null( nGHK ) && is.null( random.seed ) ) {

      result <- NextMethod( "logLik", object )

   } else {
      formula <- eval( object$call$formula )
      if( is.null( coef ) ) {
         coef <- coef( object )
      }
      if( is.null( data ) ) {
         data <- eval( object$call$data )
      }
      if( is.null( algorithm ) ) {
         algorithm <- eval( object$call$algorithm )
      }
      if( is.null( nGHK ) ) {
         nGHK <- eval( object$call$nGHK )
      }
      if( is.null( random.seed ) ) {
         random.seed <- eval( object$call$random.seed )
      }

      args <- list( formula = formula, coef = coef, data = data,
         algorithm = algorithm, nGHK = nGHK, random.seed = random.seed )
      args <- args[ !unlist( lapply( args, is.null ) ) ]

      result <- sum( do.call( mvProbitLogLik, args ) )
   }

   attr( result, "df" ) <- sum( activePar( object ) )

   class( result ) <- "logLik"

   return( result )
}
