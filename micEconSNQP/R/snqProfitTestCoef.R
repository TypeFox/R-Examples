snqProfitTestCoef <- function( nNetput, nFix, coef, form = 0,
   coefNames = c( "alpha", "beta", "gamma", "delta" ) ) {

   if( !is.list( coef ) ) {
      stop( "argument 'coef' must be a list containing the coefficients" )
   }
   if( "alpha" %in% coefNames ) {
      if( length( coef$alpha ) != nNetput ) {
         stop( "coef$alpha must have as many elements as argument 'priceNames'" )
      }
   }
   if( "beta" %in% coefNames ) {
      if( !is.matrix( coef$beta ) ) {
         stop( "coef$beta must be a matrix" )
      }
      if( nrow( coef$beta ) != ncol( coef$beta ) ) {
         stop( "coef$beta must be a _symmetric_ matrix" )
      }
      if( nrow( coef$beta ) != nNetput ) {
         stop( "coef$beta must have as many rows as argument 'priceNames' has elements" )
      }
   }
   if( "delta" %in% coefNames && nFix > 0 ) {
      if( !is.matrix( coef$delta ) ) {
         stop( "coef$delta must be a matrix" )
      }
      if( nrow( coef$delta ) != nNetput ) {
         stop( "coef$delta must have as many rows as argument 'priceNames'",
            " has elements" )
      }
      if( ncol( coef$delta ) != nFix ) {
         stop( "coef$delta must have as many columns as argument 'fixNames'",
            " has elements" )
      }
   }
   if( "gamma" %in% coefNames && nFix > 0 ) {
      if( form == 0 ) {
         if( !is.matrix( coef$gamma ) ) {
            stop( "coef$gamma must be a matrix" )
         }
         if( nrow( coef$gamma ) != ncol( coef$gamma ) ) {
            stop( "coef$gamma must be a _symmetric_ matrix" )
         }
         if( nrow( coef$gamma ) != nFix ) {
            stop( "coef$gamma must have as many rows as argument 'fixNames'",
               " has elements" )
         }
      } else if( form == 1 ) {
         if( length( dim( coef$gamma ) ) != 3 ) {
            stop( "coef$gamma must be 3-dimensional (if form == 1)" )
         }
         if( dim( coef$gamma )[ 1 ] != nNetput ) {
            stop( "the first dimension of coef$gamma must be equal to",
               " the number of elements of argument 'priceNames'" )
         }
         if( dim( coef$gamma )[ 2 ] != nFix ) {
            stop( "the second dimension of coef$gamma must be equal to",
               " the number of elements of argument 'fixNames'" )
         }
         if( dim( coef$gamma )[ 3 ] != nFix ) {
            stop( "the third dimension of coef$gamma must be equal to",
               " the number of elements of argument 'fixNames'" )
         }
      } else {
         stop( "argument 'form' must be either 0 or 1" )
      }
   }
}
