frontier <- function(
      yName, xNames = NULL, zNames = NULL, data,
      zIntercept = FALSE,
      ... ) {

   # check names of variables
   checkNames( c( yName, xNames ), names( data ) )
   if( !is.null( zNames ) ) {
      if( !is.na( zNames[1] ) ) {
         checkNames( c( zNames ), names( data ) )
      }
   }
   if( any( c( "id", "t" ) %in% c( yName, xNames, zNames ) ) ) {
      stop( "variables in arguments 'yName', 'xNames', and 'zNames'",
         " must not have names 'id' or 't'" )
   }

   # zIntercept (mu)
   if( !is.logical( zIntercept ) ) {
      stop( "argument 'zIntercept' must be logical" )
   }
   if( zIntercept &&  is.null( zNames ) ) {
      warning( "argument 'zIntercept' is ignored in",
         " Efficiency Components Frontiers (ECF)" )
   }

   # formula for the SFA
   sfaFormula <- paste( yName, "~",
      paste( c( "1", xNames ), collapse = " + " ) )

   # formula for efficiency effects
   if( is.null( zNames ) ) {
      sfaFormula <- as.formula( sfaFormula )
   } else {
      if( is.na( zNames[1] ) ) {
         if( zIntercept ) {
            effFormula <- "1"
         } else {
            effFormula <- "- 1"
         }
      } else {
         if( zIntercept ) {
            effFormula <- paste( zNames, collapse = " + " )
         } else {
            effFormula <- paste( paste( zNames, collapse = " + " ), "- 1" )
         }
      }
      sfaFormula <- as.formula( paste( sfaFormula, "|", effFormula ) )
   }

   returnObj <- sfa(
      formula = sfaFormula,
      data = data,
      ... )

   returnObj$zIntercept <- zIntercept
   returnObj$call <- match.call()
   return( returnObj )
}
