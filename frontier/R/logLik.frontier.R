logLik.frontier <- function( object, which = "mle", newParam = NULL, ... ) {

   if( is.null( newParam ) ) {
      if( tolower( which ) == "ols" ) {
         result <- object$olsLogl
      } else if( tolower( which ) == "grid" ) {
         result <- object$gridLogl
      } else if( tolower( which ) == "start" ) {
         result <- object$startLogl
      } else if( tolower( which ) == "mle" ) {
         result <- object$mleLogl
      } else {
         stop( "argument 'which' must be either 'ols', 'grid', 'start', or 'mle'" )
      }
      if( is.null( result ) ) {
         result <- NA
      }
   } else {
      if( ! is.vector( newParam ) || ! is.numeric( newParam ) ||
            length( newParam ) != length( coef( object ) ) ) {
         stop( "argument 'newParam' must be a numeric vector of length ",
            length( coef( object ) ) )
      }
      if( tolower( which ) != "mle" ) {
         warning( "argument 'which' has been ignored" )
      }
      if( "frontierQuad" %in% class( object ) ) {
         result <- logLik( frontierQuad(
            yName = eval( object$call$yName ),
            xNames = eval( object$call$xNames ),
            shifterNames = eval( object$call$shifterNames ),
            zNames = eval( object$call$zNames ),
            data = eval( object$call$data ),
            ineffDecrease = object$ineffDecrease,
            truncNorm = object$truncNorm,
            zIntercept = object$zIntercept,
            timeEffect = object$timeEffect,
            startVal = newParam,
            maxit = 0 ), which = "start" )
      } else {
         if( object$call[1] == "frontier()" ) {
            result <- logLik( frontier(
               yName = eval( object$call$yName ),
               xNames = eval( object$call$xNames ),
               zNames = eval( object$call$zNames ),
               data = eval( object$call$data ),
               ineffDecrease = object$ineffDecrease,
               truncNorm = object$truncNorm,
               zIntercept = object$zIntercept,
               timeEffect = object$timeEffect,
               startVal = newParam,
               maxit = 0 ), which = "start" )
         } else if( object$call[1] == "sfa()" ) {
            result <- logLik( sfa(
               formula = eval( object$call$formula ),
               data = eval( object$call$data ),
               ineffDecrease = object$ineffDecrease,
               truncNorm = object$truncNorm,
               timeEffect = object$timeEffect,
               startVal = newParam,
               maxit = 0 ), which = "start" )
         } else {
            stop( "unknown function '", object$call[1], "' in element 'call'" )
         }
      }
   }

   attributes( result )$nobs <- object$nob
   attributes( result )$df <- length( coef( object, which = which ) )
   class( result ) <- "logLik"

   return( result )
}