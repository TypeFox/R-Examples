cobbDouglasOpt <- function( pyName, pxNames, data, coef,
      zNames = NULL, zCoef = NULL, xNames = NULL, dataLogged = FALSE ) {

   checkNames( c( pyName, pxNames, zNames ), names( data ) )

   if( length( pyName ) != 1 | ! is.character( pyName ) ) {
      stop( "argument 'pyName' must be a single character string" )
   }

   nExog <- length( pxNames )

   if( nExog + 1 != length( coef ) ) {
      stop( "a Cobb-Douglas function with ", nExog,
         " variable exogenous variables ",
         " with corresponding prices in argument 'pxNames'",
         " must have exactly ", length( pxNames ) + 1, " coefficients",
         " in argument 'coef'" )
   }

   coefNames <- paste( "a", c( 0:nExog ), sep = "_" )
   if( is.null( names( coef ) ) ) {
      names( coef ) <- coefNames
   } else {
      coefMissing <- !( coefNames %in% names( coef ) )
      if( any( coefMissing ) ) {
         stop( "coefficient(s) ",
            paste( coefNames[ coefMissing ], collapse = ", " ),
            " are missing in argument 'coef'" )
      }
      rm( coefMissing )
   }
   rm( coefNames )

   if( length( zNames ) != length( zCoef ) ) {
      stop( "the number of fixed exogenous variables",
         " in argument 'zNames' (", length( zNames ), ")",
         " must be equal to the number of corresponding coefficients",
         " in argument 'zCoef' (", length( zCoef ), ")" )
   }

   if( !is.null( zCoef ) ) {
      zCoefNames <- paste( "d", c( 1:length( zNames ) ), sep = "_" )
      if( is.null( names( zCoef ) ) ) {
         names( zCoef ) <- zCoefNames
      } else {
         coefMissing <- !( zCoefNames %in% names( zCoef ) )
         if( any( coefMissing ) ) {
            stop( "coefficient(s) ",
               paste( coefNames[ coefMissing ], collapse = ", " ),
               " are missing in argument 'zCoef'" )
         }
         rm( coefMissing )
      }
      rm( zCoefNames )
   }

   if( is.null( xNames ) ) {
      xNames <- paste( "x", pxNames, sep = "_" )
   } else {
      if( length( xNames ) != length( pxNames ) ) {
         stop( "argument 'xNames' must have the same length as argument",
            " 'pxNames'" )
      }
   }

   if( dataLogged ) {
      logData <- data
   } else {
      logData <- logDataSet( data = data,
         varNames = c( pyName, pxNames, zNames ) )
   }

   result <- data.frame( obsNum = c( 1:nrow( logData ) ) )

   logConst <- coef[ "a_0" ]
   for( i in seq( along = zNames ) ) {
      logConst <- logConst +
         zCoef[ paste( "d", i, sep = "_" ) ] * logData[[ zNames[ i ] ]]
   }
   oneMinusSumAi <- 1 - sum( coef[ names( coef ) != "a_0" ] )

   for( i in seq( along = pxNames ) ) {
      result[[ xNames[ i ] ]] <- ( 1 / oneMinusSumAi ) * logConst
      for( j in seq( along = pxNames ) ) {
         result[[ xNames[ i ] ]] <- result[[ xNames[ i ] ]] +
            ( ( ifelse( j == i, oneMinusSumAi, 0 ) +
            coef[ paste( "a", j, sep = "_" ) ] ) / oneMinusSumAi ) *
            ( log( coef[ paste( "a", j, sep = "_" ) ] ) +
            logData[[ pyName ]] - logData[[ pxNames[ j ] ]] ) 
      }
      if( !dataLogged ) {
         result[[ xNames[ i ] ]] <- exp( result[[ xNames[ i ] ]] )
      }
   }

   if( !( "obsNum" %in% xNames ) ) {
      result$obsNum <- NULL
   }

   rownames( result ) <- rownames( data )

   return( result )
}
