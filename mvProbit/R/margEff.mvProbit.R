margEff.mvProbit <- function( object, data = eval( object$call$data ),
   cond = FALSE, othDepVar = NULL, dummyVars = object$dummyVars,
   atMean = FALSE, calcVCov = FALSE,... ) {

   # checking argument 'data'
   if( !is.data.frame( data ) ) {
      stop( "argument 'data' must be a 'data.frame'" )
   }

   # checking argument 'cond'
   if( length( cond ) != 1 ) {
      stop( "argument 'cond' must be a single logical value" )
   } else if( !is.logical( cond ) ) {
      stop( "argument 'cond' must be a logical value" )
   }

   # checking argument 'atMean'
   if( length( atMean ) != 1 ) {
      stop( "argument 'atMean' must be a single logical value" )
   } else if( !is.logical( atMean ) ) {
      stop( "argument 'atMean' must be a logical value" )
   }

   # checking argument 'calcVCov'
   if( length( calcVCov ) != 1 ) {
      stop( "argument 'calcVCov' must be a single logical value" )
   } else if( !is.logical( calcVCov ) ) {
      stop( "argument 'calcVCov' must be a logical value" )
   }

   # computing mean values of variables if requested by the user
   if( atMean ) {
      data <- as.data.frame( t( colMeans( data ) ) )
   }

   # extract the model formula
   formula <- eval( object$call$formula )

   # remove (unused) dependent variables for unconditional marginal effects
   if( !cond ) {
      formula <- formula[ - 2 ]
      if( !is.null( othDepVar ) ) {
         warning( "argument 'othDepVar' is ignored when calculating",
            " marginal effects on unconditional expectations" )
         othDepVar <- NULL
      }
   }

   # manipulate dependent variables is requested by the user
   if( !is.null( othDepVar ) ) {
      depNames <- all.vars( formula[ -3 ] )
      nDep <- length( depNames )
      nObs <- nrow( data )
      if( !is.vector( othDepVar ) || !length( othDepVar ) %in% c( 1, nDep ) ) {
         stop( "argument 'othDepVar' must be 'NULL' or a vector",
            " of length 1 or ", nDep )
      } else if( ! all( othDepVar %in% c( 0, 1, TRUE, FALSE ) ) ) {
         stop( "argument 'othDepVar' must be 'NULL' or a vector",
            " of zeros and ones of length 1 or ", nDep )
      }
      # add dependent variables that are not in the data set
      for( i in 1:nDep ) {
         if( ! depNames[ i ] %in% names( data ) ) {
            data[[ depNames[ i ] ]] <- NA
         }
      }
      for( i in 1:nObs ) {
         data[ i, depNames ] <- othDepVar
      }
   }

   if( calcVCov ) {
      vcov <- vcov( object )
   } else {
      vcov <- NULL
   }

   result <- mvProbitMargEff( formula = formula, 
      coef = coef( object ), vcov = vcov, data = data, cond = cond, 
      dummyVars = dummyVars, ... )

   return( result )
}

