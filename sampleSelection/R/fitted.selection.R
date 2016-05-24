fitted.selection <- function( object, part = "outcome", ... ) {

   if( !( part %in% c( "outcome", "selection" ) ) ) {
      stop( "argument 'part' must be either 'outcome' or 'selection'" )
   }

   # 2-step estimation
   if( object$method == "2step" ) {
      if( part == "selection" ) {
         result <- fitted( object$probit, ... )
      } else if( part == "outcome" ) {
         response <- model.frame( object$probit )[ , 1 ]
         result <- rep( NA, length( response ) )
         if( object$tobitType == 2 ) {
            result[ response == 1 ] <- fitted( object$lm, ... )
         } else if( object$tobitType == 5 ) {
            result[ response == 0 ] <- fitted( object$lm1, ... )
            result[ response == 1 ] <- fitted( object$lm2, ... )
         } else {
            stop( "unknown tobit type '",  object$tobitType,
               "' in object$tobitType" )
         }
         names( result ) <- row.names( model.frame( object$probit ) )
      } else {
         stop( "argument 'part' must be either 'outcome' or 'selection'" )
      }
   # maximum likelihood estimation
   } else if( object$method == "ml" ) {
      if( part == "selection" ) {
         modelMatrix <- model.matrix( object, part = "selection" )
         result <- drop( pnorm(modelMatrix %*% coef( object )[ object$param$index$betaS ] ) )
      } else if( part == "outcome" ) {
         if( object$tobitType == 2 ) {
            modelMatrix <- model.matrix( object )
            result <- drop(modelMatrix %*% coef( object )[ object$param$index$betaO ] )
         } else if( object$tobitType == 5 ) {
            mm <- model.matrix( object )
            response <- model.frame( object )[ , 1 ]
            result <- rep( NA, length( response ) )
            fitted1 <- mm[[ 1 ]] %*%
               coef( object )[ object$param$index$betaO1 ]
            fitted2 <- mm[[ 2 ]] %*%
               coef( object )[ object$param$index$betaO2 ]
            result[ response == 0 ] <- fitted1[ response == 0 ]
            result[ response == 1 ] <- fitted2[ response == 1 ]
         } else {
            stop( "unknown tobit type '",  object$tobitType,
               "' in object$tobitType" )
         } 
         if( object$binaryOutcome ) {
            result <- pnorm( result )
         }
      } else {
         stop( "argument 'part' must be either 'outcome' or 'selection'" )
      }
   }
   return( result )
}
