## return model matrix of the entire system
model.matrix.systemfit <- function( object, which = "x", ... ){
   result <- matrix( NA, 0, 0 )
   mmRowNames <- NULL
   mmColNames <- NULL
   for( i in 1:length( object$eq ) ) {
      mmi <- model.matrix( object$eq[[ i ]], which = which ) 
      result <- rbind(
         cbind( result, matrix( 0, nrow( result ), ncol( mmi ) ) ),
         cbind( matrix( 0, nrow( mmi ), ncol( result ) ), mmi ) )
      mmRowNames <- c( mmRowNames,
         paste( object$eq[[ i ]]$eqnLabel, "_", rownames( mmi ), sep = "" ) )
      for( j in 1:ncol( mmi ) ){
         cName <- colnames( mmi )[ j ]
         if( object$panelLike && cName != "(Intercept)" ){
            mmColNames <- c( mmColNames, cName )
         } else {
            mmColNames <- c( mmColNames,
               paste( object$eq[[ i ]]$eqnLabel, "_", cName, sep = "" ) )
         }
      }
   }
   rownames( result ) <- mmRowNames
   colnames( result ) <- mmColNames
   return( result )
}

## return model matrix of a single equation
model.matrix.systemfit.equation <- function( object, which = "x", ... ){
   if( ! which %in% c( "x", "xHat", "z" ) ) {
      stop( "argument 'which' must be either \"x\", \"xHat\", or \"z\"" )
   } else if( which %in% c( "xHat", "z" ) && is.null( object$termsInst ) ) {
      stop( "argument 'which' can only be set to \"xHat\" or \"z\" if instruments were used" )
   }
   if( which == "xHat" ) {
      xMat <- model.matrix( object, which = "x" )
      zMat <- model.matrix( object, which = "z" )
      res <- residuals( object )
      if( sum( !is.na( res ) ) != nrow( xMat ) ) {
         stop( "internal error: number of non-NA residuals not equal to",
            " number of observations in xMat. Please contact the maintainer" )
      } else if( nrow( xMat ) != nrow( zMat) ) {
         stop( "internal error: number of observations in xMat is not equal to",
            " number of observations in zMat. Please contact the maintainer" )
      }
      result <- zMat %*% solve( crossprod( zMat ), crossprod( zMat, xMat ) )
   } else {
      if( !is.null( object[[ which ]] ) ) {
         result <- object[[ which ]]
      } else if( !is.null( model.frame( object ) ) ) {
         if( which == "x" ) {
            result <- model.matrix( object$terms, data = model.frame( object ) )
         } else {
            result <- model.matrix( object$termsIns, data = object$modelInst )
         }
         attrAssign <- attributes( result )$assign
         result <- result[ !is.na( residuals( object ) ), , drop = FALSE ]
         attributes( result )$assign <- attrAssign
      } else {
         if( which == "x" ) {
            stop( "returning model matrix not possible. Please re-estimate",
               " the system with either control variable",
               " 'x' or 'model' set to TRUE" )
         } else {
            stop( "returning matrix of instruments not possible. Please re-estimate",
               " the system with either control variable",
               " 'z' or 'model' set to TRUE" )
         }
      }
   }

   return( result )
}
