## Likelihood Ratio Test
lrtest.systemfit <- function( object, ... ) {

   thisCall <- match.call()

   if( class( object ) != "systemfit" ){
      stop( "argument 'object' must be of class 'systemfit'" )
   }
   object$lrtest.systemfit.name <- deparse( substitute( object ) )
   objectList <- list( ... )
   if( length( objectList ) < 1 ){
      stop( "at least one further argument ('...') must be provided" )
   }
   if( !all( lapply( objectList, class ) == "systemfit" ) ){
      stop( "all further arguments ('...') must be of class 'systemfit'" )
   }
   dotsNames <- as.list( thisCall )[ -1 ]
   dotsNames$object <- NULL
   for( i in 1:length( objectList ) ){
      objectList[[ i ]]$lrtest.systemfit.name <- deparse( dotsNames[[ i ]] )
   }
   extractName <- function( object ){
      return( object$lrtest.systemfit.name )
   }

   result <- do.call( lrtest.default,
      c( list( object = object ), objectList, list( name = extractName ) ) )

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
