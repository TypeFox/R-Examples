prepareFixed <- function( start, activePar, fixed ) {
   nParam <- length( start )
   ## establish the active parameters.
   if(!is.null(fixed)) {
      if(!is.null(activePar)) {
         if(!all(activePar)) {
            warning("Both 'activePar' and 'fixed' specified.  'activePar' ignored")
         }
      }
      if( is.logical( fixed ) ) {
         if( length ( fixed ) != length( start ) || !is.null( dim( fixed ) ) ) {
            stop( "if fixed parameters are specified using logical values,",
               " argument 'fixed' must be a logical vector",
               " with one element for each parameter",
               " (number of elements in argument 'start')" )
         }
         activePar <- !fixed
      } else if( is.numeric( fixed ) ) {
         if( length ( fixed ) >= length( start ) || !is.null( dim( fixed ) ) ) {
            stop( "if fixed parameters are specified using their positions,",
               " argument 'fixed' must be a numerical vector",
               " with less elements than the number of parameters",
               " (number of elements in argument 'start'" )
         } else if( min( fixed ) < 1 || max(fixed ) > length( start ) ) {
            stop( "if fixed parameters are specified using their positions,",
               " argument 'fixed' must have values between 1 and",
               " the total number of parameter",
               " (number of elements in argument 'start'" )
         }
         activePar <- ! c( 1:length( start ) ) %in% fixed
      } else if( is.character( fixed ) ) {
         if( length ( fixed ) >= length( start ) || !is.null( dim( fixed ) ) ) {
            stop( "if fixed parameters are specified using their names,",
               " argument 'fixed' must be a vector of character strings",
               " with less elements than the number of parameters",
               " (number of elements in argument 'start'" )
         } else if( is.null( names( start ) ) ) {
            stop( "if fixed parameters are specified using their names,",
               " parameter names have to be specified in argument 'start'" )
         } else if( any( ! names( fixed ) %in% names( start ) ) ) {
            stop( "if fixed parameters are specified using their names,",
               " all parameter names specified in argument 'fixed'",
               " must be specified in argument 'start'" )
         }
         activePar <- ! names( start ) %in% fixed
      } else {
         stop( "argument 'fixed' must be either a logical vector,",
            " a numeric vector, or a vector of character strings" )
      }
   } else {
      if( is.null( activePar ) ) {
         activePar <- rep( TRUE, length( start ) )
      } else if(is.numeric(activePar)) {
         a <- rep(FALSE, nParam)
         a[activePar] <- TRUE
         activePar <- a
      }
   }
   names( activePar ) <- names( start )
   if( all( !activePar ) ){
      stop( "At least one parameter must not be fixed",
         " using argument 'fixed'" )
   }
   return( !activePar )
}