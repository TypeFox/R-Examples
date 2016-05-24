cesEstStart <- function( yName, xNames, tName, data, vrs,
      method, start, rho1, rho2, rho, nParam, nested = FALSE,
      multErr ) {

   withTime <- !is.null( tName )

   # number of explanatory variables
   nExog <- length( xNames )

   # start values
   if( method %in% c( "Kmenta", "DE" ) ) {
      if( !is.null( start ) ) {
         warning( "ignoring starting values because they are not required",
            " for method '", method, "'" )
         start <- NULL
      }
   } else {
      if( is.null( start ) ) {
         rhoStart <- ifelse( is.null( rho ), 0.25, rho )
         rho1Start <- ifelse( is.null( rho1 ), 0.25, rho1 )
         rho2Start <- ifelse( is.null( rho2 ), 0.25, rho2 )
         if( nested && nExog == 3 ) {
            start <- c( 1, 0.5, 0.5, rho1Start, rhoStart )
         } else if( nested && nExog == 4 ) {
            start <- c( 1, 0.5, 0.5, 0.5, rho1Start, rho2Start, rhoStart )
         } else if( !nested && nExog == 2 ) {
            start <- c( 1, 0.5, rhoStart )
         } else {
            stop( "cannot create starting values for a",
               ifelse( nested, " nested", " non-nested" ), " CES function",
               " with ", nExog, " explanatory variables" )
         }
         if( withTime ) {
            start <- c( start[ 1 ], 0.015, start[ -1 ] )
         }
         if( vrs ) {
            start <- c( start, 1 )
         }
         yTemp <- cesCalc( xNames = xNames, data = data, coef = start,
            tName = tName, nested = nested )
         if( multErr ) {
            start[ 1 ] <- 
               mean( log( data[[ yName ]] ), na.rm = TRUE ) -
               mean( log( yTemp ), na.rm = TRUE )
         } else {
            start[ 1 ] <- 
               mean( data[[ yName ]], na.rm = TRUE ) /
               mean( yTemp, na.rm = TRUE )
         }
         start[ 1 ] <- max( start[ 1 ], 1e-8 )
         if( !is.null( rho ) ) {
            start <- start[ -( withTime + ifelse( nested, 2 * nExog - 1, 3 ) ) ]
         }
         if( !is.null( rho2 ) && nested && nExog == 4 ) {
            start <- start[ -( 6 + withTime ) ]
         }
         if( !is.null( rho1 ) && nested ) {
            start <- start[ -( 5 - ( nExog == 3 ) + withTime ) ]
         }
      }
      if( length( start ) != nParam ) {
         stop( "wrong number of starting values:",
            " you provided ", length( start ), " values",
            " but the model has ", nParam, " parameters" )
      }
      names( start ) <- cesCoefNames( nExog = length( xNames ), vrs = vrs,
         returnRho = is.null( rho ), returnRho1 = is.null( rho1 ), 
         returnRho2 = is.null( rho2 ), nested = nested, withTime = withTime )
      # checking starting values
      if( any( is.infinite( start ) ) ) {
         stop( "all starting values must be finite" )
      }
      # checking gamma
      if( nested && nExog == 3 ) {
         if( start[ "gamma" ] <= 0 ) {
            stop( "the starting value for 'gamma' must be positive" )
         }
      } else {
         if( start[ "gamma" ] <= 0 ) {
            stop( "the starting value for 'gamma' must be positive" )
         }
      }
      # checking delta
      if( nested ) {
         if( start[ "delta_1" ] < 0 || start[ "delta_1" ] > 1 ) {
            stop( "the starting value for 'delta_1' must be between 0 and 1" )
         }
         if( nExog == 4 ) {
            if( start[ "delta_2" ] < 0 || start[ "delta_2" ] > 1 ) {
               stop( "the starting value for 'delta_2' must be between 0 and 1" )
            }
         }
      }
      if( start[ "delta" ] < 0 || start[ "delta" ] > 1 ) {
            stop( "the starting value for 'delta' must be between 0 and 1" )
      }
      # checking rho
      if( is.null( rho ) ) {
         if( start[ "rho" ] < -1 ) {
            stop( "the starting value for 'rho' must be -1 or larger" )
         }
      }
      if( is.null( rho1 ) && nested ) {
         if( start[ "rho_1" ] < -1 ) {
            stop( "the starting value for 'rho_1' must be -1 or larger" )
         }
      }
      if( is.null( rho2 ) && nested && nExog == 4 ) {
         if( start[ "rho_2" ] < -1 ) {
            stop( "the starting value for 'rho_2' must be -1 or larger" )
         }
      }
      # checking nu
      if( vrs ) {
         if( start[ "nu" ] <= 0 ) {
            stop( "the starting value for 'nu' must be positive" )
         }
      }
   }

   return( start )
}
