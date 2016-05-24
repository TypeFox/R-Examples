.micEconCoefOrder <- function( coefNames ) {

   nCoef <- length( coefNames )

   coefChar <- rep( NA, nCoef )
   coefNum1 <- rep( NA, nCoef )
   coefNum2 <- rep( NA, nCoef )

   parts <- strsplit( coefNames, "_" )

   for( i in 1:nCoef ) {
      coefChar[ i ] <- parts[[ i ]][ 1 ]
      coefNum1[ i ] <- as.integer( parts[[ i ]][ 2 ] )
      if( length( parts[[ i ]] ) > 2 ) {
         if( regexpr( "^[0-9]*$", parts[[ i ]][ 3 ] ) != -1 ) {
            coefNum2[ i ] <- as.integer( parts[[ i ]][ 3 ] )
         }
      }
   }

   result <- order( coefChar, as.integer( coefNum1 ), as.integer( coefNum2 ) )

   return( result )
}