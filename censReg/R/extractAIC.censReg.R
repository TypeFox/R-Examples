extractAIC.censReg <- function( fit, scale = 0, k = 2, ... ) {
   # copied from stats:::extractAIC.glm and slightly modified thereafter
   n <- nObs( fit )
   edf <- n - df.residual( fit )
   aic <- AIC( fit )
   return( c( edf, aic + ( k - 2 ) * edf ) )
}

