cesCheckRhoApprox <- function( rhoApprox, withY, withDeriv ) {


   if( !is.vector( rhoApprox ) || !is.numeric( rhoApprox ) ) {
      stop( "argument 'rhoApprox' must be a numeric vector" )
   } 

   if( is.na( withY ) && is.na( withDeriv ) ) {
      stop( "internal error: arguments 'withY' and 'withDeriv'",
         " must not be 'NA' at the same time" )
   } else if( is.na( withY ) ) {
      if( length( rhoApprox ) == 5 ) {
         rhoApprox <- rhoApprox[ -1 ]
      }
      withY <- FALSE
      if( length( rhoApprox ) != 4 ) {
         stop( "argument 'rhoApprox' must be a vector of length",
            " 4 or 5 but it has a length of ",
            length( rhoApprox ) )
      }
   } else if( is.na( withDeriv ) ) {
      if( length( rhoApprox ) == 5 ) {
         rhoApprox <- rhoApprox[ 1 ]
      }
      withDeriv <- FALSE
      if( length( rhoApprox ) != 1 ) {
         stop( "argument 'rhoApprox' must be a vector of length",
            " 1 or 5 but it has a length of ",
            length( rhoApprox ) )
      }
   } else {
      if( length( rhoApprox ) != ( withY + 4 * withDeriv ) ) {
         stop( "argument 'rhoApprox' must be a vector of length ",
            withY + 4 * withDeriv, " but it has a length of ",
            length( rhoApprox ) )
      }
   }

   elemNames <- NULL
   if( withY ) {
      elemNames <- c( elemNames, "y" )
   }
   if( withDeriv ) {
      elemNames <- c( elemNames, c( "gamma", "delta", "rho", "nu" ) )
   }

   if( !is.null( names( rhoApprox ) ) ) {
      if( any( names( rhoApprox ) != elemNames ) ) {
         warning( "ignoring names of elements in the vector",
            " that was provided as argument 'rhoApprox'" )
      }
   }
   names( rhoApprox ) <- elemNames

   return( rhoApprox )
}

