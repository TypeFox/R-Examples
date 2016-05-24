print.summary.cesEst <- function( x, ela = TRUE,
      digits = max( 3, getOption( "digits" ) - 3 ), ... ) {

   cat( "Estimated CES function with " )
   if( "nu" %in% rownames( coef( x ) ) ){
      cat( "variable " )
   } else {
      cat( "constant " )
   }
   cat( "returns to scale\n\n" )

   cat( "Call:\n" )
   print( x$call )
   cat( "\n" )

   cat( "Estimation by " )
   if( x$method == "Kmenta" ) {
      cat( "the linear Kmenta approximation\n" )
      cat( "Test of the null hypothesis that the restrictions of the Translog\n",
         "function required by the Kmenta approximation are true:\n",
         "P-value = ", x$testKmenta[ 2, "Pr(>F)" ], "\n", sep = "" )
   } else {
      cat( "non-linear least-squares using the '", x$method, "' optimizer\n",
         sep = "" )
      if( !is.null( x$allRhoSum ) ) {
         gridCoef <- c( "rho1", "rho2", "rho" ) 
         gridCoef <- gridCoef[ gridCoef %in% names( x$allRhoSum ) ]
         cat( "and a ",
            c( "one", "two", "three" )[ length( gridCoef ) ],
            "-dimensional grid search for coefficient",
            ifelse( length( gridCoef ) > 1, "s ", " " ),
            paste( 
               paste( "'",sub( "([12])$", "_\\1", gridCoef ), "'", sep = "" ), 
               collapse = ", " ),
            "\n", sep = "" )
      }
      cat( "assuming",
         ifelse( x$multErr, "a multiplicative", "an additive" ),
         "error term\n" )
      if( is.null( x$allRhoSum ) ) {
         if( !is.null( x[[ "rho1" ]] ) ) {
            cat( "Coefficient 'rho_1' was fixed at", x$rho1, "\n" )
         }
         if( !is.null( x[[ "rho2" ]] ) ) {
            cat( "Coefficient 'rho_2' was fixed at", x$rho2, "\n" )
         }
         if( !is.null( x[[ "rho" ]] ) ) {
            cat( "Coefficient 'rho' was fixed at", x$rho, "\n" )
         }
      }
      if( !is.null( x$convergence ) ) {
         cat( "Convergence ", ifelse( x$convergence, "", "NOT " ),
            "achieved after ", sep = "" )
         if( length( x$iter ) == 1 ) {
            cat( x$iter, "iterations\n" )
         } else {
            cat( paste( x$iter, names( x$iter ), collapse = " and " ),
               "calls\n" )
         }
      }
      if( !is.null( x$message ) ) {
         cat( "Message:", x$message, "\n" )
      }
   }
   cat( "\n" )

   cat( "Coefficients:\n" )
   printCoefmat( coef( x ), digits = digits )
   cat( "\n" )

   cat( "Residual standard error:", x$sigma, "\n" )
   cat( "Multiple R-squared:", x$r.squared, "\n" )
   cat( "\n" )

   if( ela && is.matrix( x$ela ) ) {
      if( nrow( x$ela ) == 1 ) {
         cat( "Elasticity of Substitution:\n" )
      } else {
         cat( "Elasticities of Substitution:\n" )
      }
      printCoefmat( x$ela, digits = digits )
      if( nrow( x$ela ) > 1 ) {
         cat( "HM = Hicks-McFadden (direct) elasticity of substitution\n" )
         cat( "AU = Allen-Uzawa (partial) elasticity of substitution\n" )
      }
      cat( "\n" )
   }

   invisible( x )
}
