.aidsCheckCoef <- function( coef, argCoef = "coef",
      variables = NULL ) {

   # checking argument 'variables' of *this* function
   if( !is.null( variables ) ){
      if( !is.list( variables ) ){
         stop( "internal error: argument 'variables' must be of class 'list'" )
      }
      for( i in 1:length( variables ) ){
         if( !is.list( variables[[ i ]] ) ){
            stop( "internal error: each element of argument 'variables'",
               " must be of class 'list'" )
         }
         if( length( variables[[ i ]] ) != 3 ){
            stop( "internal error: each element of argument 'variables'",
               " must have 3 elements" )
         }
         if( !( is.numeric( variables[[ i ]][[ 1 ]] ) || 
               is.na( variables[[ i ]][[ 1 ]] ) ) || 
               length( variables[[ i ]][[ 1 ]] ) != 1 ){
            stop( "internal error: 'variables[[ ", i , " ]][[ 1 ]]' must be",
               " a numeric scalar" )
         }
         if( !is.character( variables[[ i ]][[ 2 ]] ) ){
            stop( "internal error: 'variables[[ ", i , " ]][[ 2 ]]' must be",
               " a character string" )
         }
         if( !( variables[[ i ]][[ 3 ]] %in% c( "goods", "shifters" ) ) ){
            stop( "internal error: 'variables[[ ", i , " ]][[ 3 ]]' must be",
               " either 'goods' or 'shifters'" )
         }
      }
   }

   ## checking coefficients
   # alpha
   if( is.null( coef$alpha ) ){
      return( paste( "'", argCoef, "$alpha' is missing", sep = "" ) )
   }
   if( !is.numeric( coef$alpha ) ){
      return( paste( "'", argCoef, "$alpha' must be numeric vector", sep = "" ) )
   }

   # beta
   if( is.null( coef$beta ) ){
      return( paste( "'", argCoef, "$beta' is missing", sep = "" ) )
   }
   if( !is.numeric( coef$beta ) ){
      return( paste( "'", argCoef, "$beta' must be numeric vector", sep = "" ) )
   }
   if( length( coef$alpha ) != length( coef$beta ) ) {
      return( paste( "'", argCoef, "$alpha' and '", argCoef, "$beta'",
         " must have the same length", sep = "" ) )
   }

   # gamma
   if( is.null( coef$gamma ) ){
      return( paste( "'", argCoef, "$gamma' is missing", sep = "" ) )
   }
   if( !is.matrix( coef$gamma ) ){
      return( paste( "'", argCoef, "$gamma' must be a matrix", sep = "" ) )
   }
   if( nrow( coef$gamma ) != ncol( coef$gamma ) ) {
      return( paste( "argument '", argCoef, "$gamma' must be a square matrix",
         sep = "" ) )
   }
   if( length( coef$alpha ) != nrow( coef$gamma ) ) {
      return( paste( "the number of rows of '", argCoef, "$gamma'",
         " must be equal to the length of '", argCoef, "$alpha'",
         sep = "" ) )
   }

   # delta
   if( !is.null( coef$delta ) ){
      if( !is.matrix( coef$delta ) ){
         return( paste( "'", argCoef, "$delta' must be a matrix", sep = "" ) )
      }
      if( length( coef$alpha ) != nrow( coef$delta ) ) {
         return( paste( "the number of rows of '", argCoef, "$delta'",
            " must be equal to the length of '", argCoef, "$alpha'",
            sep = "" ) )
      }
   }

   # checking variables
   if( !is.null( variables ) ){
      for( i in 1:length( variables ) ){
         if( variables[[ i ]][[ 3 ]] == "goods" ){ 
            if( variables[[ i ]][[ 1 ]] != length( coef$alpha ) && 
                  !is.na( variables[[ i ]][[ 1 ]] ) ) {
               return( paste( "'", argCoef, "$alpha' and '", variables[[ i ]][[ 2 ]],
                  "' must have the same length", sep = "" ) )
            }
         } else if( variables[[ i ]][[ 3 ]] == "shifters" ){
            if( variables[[ i ]][[ 1 ]] != ncol( coef$delta ) && 
                  !is.na( variables[[ i ]][[ 1 ]] ) ) {
               return( paste( "the number of columns of '", argCoef, "$delta'",
                  " must be equal to the length of '", variables[[ i ]][[ 2 ]], "'", 
                  sep = "" ) )
            }
         } else {
            stop( "internal error: 'variables[[ ", i , " ]][[ 3 ]]' must be",
               " either 'goods' or 'shifters'" )
         }
      }
   }

   return( NULL )
}