## Likelihood Ratio Test
lrtest.aidsEst <- function( object, ... ) {

   thisCall <- match.call()
   dotsList <- list( ... )
   dotsNames <- as.list( thisCall )[ -1 ]
   dotsNames$object <- NULL

   if( length( dotsList ) < 1 ){
      stop( "at least two arguments (objects of class 'aidsEst')",
         " are required" )
   }
   if( class( object ) != "aidsEst" ||
         !all( lapply( dotsList, class ) == "aidsEst" ) ){
      stop( "all arguments must be of class 'aidsEst'" )
   }

   createLabel <- function( aidsEstObject, objectName ){
      if( aidsEstObject$method == "LA" ) {
         result <- "LA-AIDS"
      } else {
         result <- "AIDS"
      }
      if( is.null( aidsEstObject$call$sym ) ) {
         aidsEstObject$call$sym <- TRUE
      }
      if( is.null( aidsEstObject$call$hom ) ) {
         aidsEstObject$call$hom <- TRUE
      }
      if( aidsEstObject$call$sym ) {
         result <- paste( result, ", symmetry and homogeneity imposed", sep = "" )
      } else if( aidsEstObject$call$hom ) {
         result <- paste( result, ", homogeneity imposed", sep = "" )
      } else {
         result <- paste( result, ", unrestricted", sep = "" )
      }
      if( !is.null( aidsEstObject$call$shifterNames ) ) {
         result <- paste( result, ", ",
            length( aidsEstObject$call$shifterNames ) - 1,
            " demand shifter(s)", sep = "" )
      }
      result <- paste( objectName, " (", result, ")", sep = "" )
      return( result )
   }

   object$lrtest.aidsEst.name <- createLabel( object,
      deparse( substitute( object ) ) )

   for( i in 1:length( dotsList ) ){
      dotsList[[ i ]]$lrtest.aidsEst.name <- createLabel( dotsList[[ i ]],
         deparse( dotsNames[[ i ]] ) )
   }
   extractName <- function( object ){
      return( object$lrtest.aidsEst.name )
   }

   result <- do.call( lrtest.default,
      c( list( object = object ), dotsList, list( name = extractName ) ) )

   for( i in 2:nrow( result ) ){
      if( ( result[ i, "#Df" ] - result[ i - 1, "#Df" ] ) *
            ( result[ i, "LogLik" ] - result[ i - 1, "LogLik" ] ) < 0 ) {
         if( result[ i, "LogLik" ] > result[ i - 1, "LogLik" ] ) {
            compareLikelihood <- "larger"
            compareDf <- "less"
         } else {
            compareLikelihood <- "smaller"
            compareDf <- "more"
         }
         warning( "model '", i, "' has a ", compareLikelihood,
            " log-likelihood value than the ", compareDf,
            " restricted model '", i - 1, "'" )
      }
   }

   return( result )
}
