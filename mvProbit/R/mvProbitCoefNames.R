mvProbitCoefNames <- function( nDep, nReg ) {
   result <- c( 
      paste( "b", rep( 1:nDep, each = nReg ), rep( 0:(nReg-1), nDep ), 
         sep = "_" ),
      matrix( paste( "R", rep( 1:nDep, each = nDep ), rep( 1:nDep, nDep ),
         sep = "_" ), nrow = nDep )[ lower.tri( diag(nDep) ) ] )
   return( result )
}