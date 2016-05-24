quadFuncEla <- function( xNames, data, coef,
   yName = NULL, shifterNames = NULL, homWeights = NULL ) {

   checkNames( c( xNames, yName ), names( data ) )

   # if 'data' is a vector, convert it to a data.frame
   data <- .micEconVectorToDataFrame( data )

   # check argument 'homWeights'
   .quadFuncCheckHomWeights( homWeights, xNames )

   # check argument 'coef'
   .quadFuncCheckCoefNames( names( coef ), nExog = length( xNames ),
      shifterNames = shifterNames, data = data, warn = FALSE )

   if( is.null( yName ) ){
      yHat <- quadFuncCalc( xNames = xNames, data = data, coef = coef, 
         shifterNames = shifterNames, homWeights = homWeights )
   } else {
      yHat <- data[[ yName ]]
   }

   result <- quadFuncDeriv( xNames = xNames, data = data, coef = coef,
      homWeights = homWeights )
   for( i in xNames ) {
      result[[ i ]] <- result[[ i ]] * data[[ i ]] / yHat
   }

   class( result ) <- c( "quadFuncEla", class( result ) )
   return( result )
}
