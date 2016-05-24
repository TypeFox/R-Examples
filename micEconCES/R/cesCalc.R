cesCalc <- function( xNames, data, coef, tName = NULL,
      nested = FALSE, rhoApprox = 5e-6 ) {

   # check number of exogenous variables
   nExog <- length( xNames )
   if( nExog < 2 ) {
      stop( "argument 'xNames' must include the names of at least 2 variables" )
   }

   # check time variable
   if( is.null( tName ) ) {
      withTime <- FALSE
   } else {
      withTime <- TRUE
      if( length( tName ) != 1 || !is.character( tName[1] ) ) {
         stop( "argument 'tName' must be a single character string" )
      }
   }

   # check 
   if( nested && ! nExog %in% c( 3, 4 ) ) {
      stop( "the nested CES function is currently implemented only for",
         " 3 and 4 inputs" )
   }

   # check names of exogenous variables
   checkNames( c( xNames, tName ), names( data ) )

   # check argument 'rhoApprox'
   if( !is.numeric( rhoApprox ) || length( rhoApprox ) != 1 ) {
      stop( "argument 'rhoApprox' must be a numeric scalar" )
   }

   # check for VRS
   if( nExog == 2 ) {
      vrs <- length( coef ) >= 4 + withTime
   } else {
      if( nested ) {
         vrs <- length( coef ) >= 2 * nExog + withTime
      } else {
         vrs <- length( coef ) - nExog  >= 3 + withTime
      }
   }

   # check number of coefficients
   if( nExog == 2 && length( coef ) != 3 + vrs + withTime ) {
      stop( "a CES function with 2 exogenous variables and",
         ifelse( vrs, " variable", " constant" ), " returns to scale",
         ifelse( withTime, " and time effect", "" ),
         " must have ", 3 + vrs + withTime, " coefficients",
         " but you provided ", length( coef ), " coefficients" )
   } else if( nExog > 2 && !nested && length( coef ) != nExog + 2 + vrs + withTime ) {
      stop( "a non-nested CES function with ", nExog, " exogenous variables and",
         ifelse( vrs, " variable", " constant" ), " returns to scale",
         ifelse( withTime, " and time effect", "" ),
         " must have ", nExog + 2 + vrs + withTime, " coefficients",
         " but you provided ", length( coef ), " coefficients" )
   } else if( nExog %in% c( 3, 4 ) && nested && length( coef ) != 
         2 * nExog - 1 + vrs + withTime ) {
      stop( "a nested CES function with ", nExog, " exogenous variables and",
         ifelse( vrs, " variable", " constant" ), " returns to scale",
         ifelse( withTime, " and time effect", "" ),
         " must have ", 2 * nExog - 1 + vrs + withTime, " coefficients",
         " but you provided ", length( coef ), " coefficients" )
   }

   # check for NAs in coefficients
   if( sum( is.na( coef ) ) > 0 ) {
      warning( "some of the coefficiencients are 'NA'" )
   }

   # names of coefficients
   coefNames <- cesCoefNames( nExog = nExog, vrs = vrs, nested = nested,
      withTime = withTime )

   # assign or check names of coefficients
   if( is.null( names( coef ) ) ) {
      names( coef ) <- coefNames
   } else {
      coefTest <- coefNames %in% names( coef )
      if( any( !coefTest ) ) {
         stop( "following coefficient name(s) is/are missing in argument",
            " 'coef': ", paste( coefNames[ !coefTest ], collapse = ", " ) )
      }
   }

   # make the case of two explanatory compatible to the case of N variables
   if( nExog == 2 ) {
      names( coef )[ names( coef ) == "delta" ] <- "delta_1"
      coef <- c( coef, delta_2 = 1 - unname( coef[ "delta_1" ] ) )
   }

   # check if the deltas sum up to one
   if( !nested ) {
      deltaCoefs <- coef[ grep( "delta\\_", names( coef ) ) ]
      if( sum( is.na( deltaCoefs ) ) == 0 &&
            abs( sum( deltaCoefs, na.rm = TRUE ) - 1 ) / sum( abs( deltaCoefs ) ) >
            .Machine$double.eps^0.5 ) {
         stop( "the sum of the delta coefficients must sum up to 1" )
      }
   }

   # make the case of constant returns to scale (CRS) compatible to the VRS case
   if( !vrs ) {
      coef <- c( coef, nu = 1 )
   }

   # calculate the endogenous variable
   if( !nested ) {
      gammaStar <- coef[ "gamma" ]
      if( withTime ) {
         gammaStar <- gammaStar * exp( coef[ "lambda" ] * data[[ tName ]] )
      } 
      if( abs( coef[ "rho" ] ) <= rhoApprox ) {
         result <- log( gammaStar )
         for( i in 1:nExog ) {
            result <- result + coef[ paste( "delta", i, sep = "_" ) ] *
               coef[ "nu" ] * log( data[[ xNames[ i ] ]] )
         }
         for( i in 1:( nExog - 1 ) ) {
            for( j in ( i + 1 ):nExog ) {
               result <- result - 0.5 * coef[ "rho" ] * coef[ "nu" ] *
                  coef[ paste( "delta", i, sep = "_" ) ] *
                  coef[ paste( "delta", j, sep = "_" ) ] *
                  ( log( data[[ xNames[ i ] ]] ) - log( data[[ xNames[ j ] ]] ) )^2
            }
         }
         result <- exp( result )
      } else {
         if( coef[ "rho" ] == 0 ) {
            result <- NaN
         } else {
            result <- 0
            for( i in 1:nExog ) {
               result <- result + coef[ paste( "delta", i, sep = "_" ) ] *
                  data[[ xNames[ i ] ]]^( -coef[ "rho" ] )
            }
            result <- result^( -coef[ "nu" ] / coef[ "rho" ] )
            result <- gammaStar * result
         }
      }
   } else if( nExog == 3 ) {   # nested CES with 3 inputs
         result <- cesInterN3( funcName = "cesCalcN3", 
            par = coef, xNames = xNames, tName = tName,
            data = data, rhoApprox = rhoApprox )
   } else {                    # nested CES with 4 inputs
      result <- cesInterN4( funcName = "cesCalcN4", 
            par = coef, xNames = xNames, tName = tName,
            data = data, rhoApprox = rhoApprox )
   }

   return( result )
}
