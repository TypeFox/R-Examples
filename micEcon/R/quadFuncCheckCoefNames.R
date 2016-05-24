.quadFuncCheckCoefNames <- function( coefNames, nExog,
      shifterNames = NULL, data = NULL, warn = TRUE ) {

   if( is.null( coefNames ) ) {
      stop( "argument 'coef' does not have names" )
   }

   # alphas
   allCoefNames <- paste( "a", 0:nExog, sep = "_" )

   # betas
   if( nExog > 0 ) {
      for( i in 1:nExog ) {
         for( j in i:nExog ) {
            allCoefNames <- c( allCoefNames, paste( "b", i, j, sep = "_" ) )
         }
      }
   }

   # deltas
   for( i in seq( along = shifterNames ) ) {
      if( is.logical( data[[ shifterNames[ i ] ]] ) ) {
         allCoefNames <- c( allCoefNames, paste( "d", i, "TRUE", sep = "_" ) )
      } else if( is.factor( data[[ shifterNames[ i ] ]] ) ) {
         noOmittedLevels <- 0
         for( j in levels( data[[ shifterNames[ i ] ]] ) ) {
            thisCoefName <- paste( "d", i, j, sep = "_" )
            if( thisCoefName %in% names( coef ) | noOmittedLevels > 0 ) {
               allCoefNames <- c( allCoefNames, thisCoefName )
            } else {
               noOmittedLevels <- noOmittedLevels + 1
            }
         }
      } else {
         allCoefNames <- c( allCoefNames, paste( "d", i, sep = "_" ) )
      }
   }

   missingNames <- allCoefNames[ ! allCoefNames %in% coefNames ]

   if( length( missingNames ) > 0 ) {
      stop( "following coefficients in argument 'coef' are missing: ",
         paste( missingNames, collapse = ", " ) )
   }

   if( warn ) {
      coefNames[ coefNames == "" ] <-
         c( 1:length( coefNames ) )[ coefNames == "" ]
      nonUsedCoefs <- coefNames[ ! coefNames %in% allCoefNames ]

      if( length( nonUsedCoefs ) > 0 ) {
         warning( "following coefficients in argument 'coef' are not used: ",
            paste( nonUsedCoefs, collapse = ", " ) )
      }
   }

   return( NULL )
}