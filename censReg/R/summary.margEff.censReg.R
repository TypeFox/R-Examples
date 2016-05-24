summary.margEff.censReg <- function( object, ... ) {
   
   vcov <- attr( object, "vcov" )
   if( is.null( vcov ) ) {
      warning( "cannot calculate standard errors, t-values, and P-values,",
         " because the marginal effects do not have an attribute 'vcov':",
         " please set attribute 'calcVCov' of the 'margEff' method to 'TRUE'" )
      stdEr <- rep( NA, length( object ) )
   } else {
      stdEr <- sqrt( diag( vcov ) )
   }
   tStat <- object / stdEr
   pVal <- 2 * ( 1 - pt( abs( tStat ), attr( object, "df.residual" ) ) )
   result <- cbind( object, stdEr, tStat, pVal )
   colnames( result ) <- 
      c( "Marg. Eff.", "Std. Error", "t value", "Pr(>|t|)" )

   class( result ) <- c( "summary.margEff.censReg", class( result ) )

   return( result )
}