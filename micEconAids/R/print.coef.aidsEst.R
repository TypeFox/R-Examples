print.coef.aidsEst <- function( x, ... ) {
   if( !is.null( x$alpha0 ) ){
      cat( "alpha_0\n" )
      cat( x$alpha0, "\n" )
   }
   cat( "alpha\n" )
   print( x$alpha )
   cat( "beta\n" )
   print( x$beta )
   cat( "gamma\n" )
   print( x$gamma )
   if( !is.null( x$delta ) ){
      cat( "delta\n" )
      print( x$delta )
   }
   invisible( x )
}
