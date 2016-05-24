model.matrix.selection <- function( object, part = "outcome",
                                   ... ) {
   ## part   'outcome' or 'selection'
   ##        find the design matrix for that submodel
   ##        
   if( !( part %in% c( "outcome", "selection" ) ) ) {
      stop( "argument 'part' must be either 'outcome' or 'selection'" )
   }

   # 2-step estimation
   if( object$method == "2step" ) {
      if( part == "selection" ) {
         result <- model.matrix( object$probit, ... )
      } else if( part == "outcome" ) {
         response <- model.frame( object$probit )[ , 1 ]
         nObs <- length( response )
         obsNames <- row.names( model.frame( object$probit ) )
         if( tobitType(object) == 2 ) {
            mm <- model.matrix( object$lm, ... )
            result <- matrix( NA, nrow = nObs, ncol = ncol( mm ) )
            result[ response == 1, ] <- mm
            attributes( result )$assign <- attributes( mm )$assign
            attributes( result )$contrasts <- attributes( mm )$contrasts
            rownames( result ) <- obsNames
            colnames( result ) <- colnames( mm )
         } else if( object$tobitType == 5 ) {
            result <- list()
            mm <- list()
            mm[[ 1 ]] <- model.matrix( object$lm1, ... )
            mm[[ 2 ]] <- model.matrix( object$lm2, ... )
            for( i in 1:2 ) {
               result[[ i ]] <- matrix( NA, nrow = nObs, ncol = ncol( mm[[ i ]] ) )
               result[[ i ]][ response == ( i - 1 ), ] <- mm[[ i ]]
               attributes( result[[ i ]] )$assign <-
                  attributes( mm[[ i ]] )$assign
               attributes( result[[ i ]] )$contrasts <-
                  attributes( mm[[ i ]] )$contrasts
               rownames( result[[ i ]] ) <- obsNames
               colnames( result[[ i ]] ) <- colnames( mm[[ i ]] )
            }
         } else {
            stop( "unknown tobit type '",  object$tobitType,
               "' in object$tobitType" )
         }
      } else {
         stop( "argument 'part' must be either 'outcome' or 'selection'" )
      }
   }
   ## maximum likelihood estimation
   else if( object$method == "ml" ) {
      if( part == "selection" ) {
         if( ! is.null( object$xs ) ) {
            result <- object$xs
         } else {
            mf <- model.frame( object )
            result <- model.matrix( object$termsS, mf)
         }
      }
      else if( part == "outcome" ) {
         response <- model.frame( object )[ , 1 ]
         nObs <- length( response )
         if( object$tobitType == 2 ) {
            if( ! is.null( object$xo ) ) {
               result <- object$xo
            } else {
               mf <- model.frame( object )
               attributes( mf )$na.action <- na.pass
               result <- model.matrix( object$termsO, mf )
               result[ response == 0, ] <- NA
            }
         }
         else if( object$tobitType == 5 ) {
            result <- list("1"=NULL, "2"=NULL)
            if( ! is.null( object$xo1 ) && ! is.null( object$xo2 ) ) {
               result[[ 1 ]] <- object$xo1
               result[[ 2 ]] <- object$xo2
            } else {
               mf <- model.frame( object )
               attributes( mf )$na.action <- na.pass
               result[[ 1 ]] <- model.matrix( object$termsO1, mf )
               result[[ 2 ]] <- model.matrix( object$termsO2, mf )
               result[[ 1 ]][ response == 1, ] <- NA
               result[[ 2 ]][ response == 0, ] <- NA
            }
         }
         else if(tobitType(object) == "treatment" ) {
            ## It is like tobit5, just one outcome, and everything is
            ## observed 
              if( ! is.null( object$xo ) ) {
                 result <- object$xo
              }
              else {
                 mf <- model.frame( object, ... )
                           # '...' passes 'data' argument
                 attributes( mf )$na.action <- na.pass
                 result <- model.matrix( object$termsO, mf)
              }
         }
      }
   }

   return( result )
}
