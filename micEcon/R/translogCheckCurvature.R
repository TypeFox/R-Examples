translogCheckCurvature <- function( xNames, data, coef, convexity = TRUE,
   quasi = FALSE, dataLogged = FALSE, ... ) {

   result <- list()

   hessian <- translogHessian( xNames = xNames, data = data, coef = coef,
      dataLogged = dataLogged, bordered = quasi )

   if( quasi ) {
      if( convexity ) {
         result$obs <- quasiconvexity( hessian, ... )
      } else {
         result$obs <- quasiconcavity( hessian, ... )
      }
   } else {
      result$obs <- semidefiniteness( hessian, positive = convexity, ... )
   }

   result$convexity <- convexity
   result$quasi     <- quasi

   class( result ) <- "translogCheckCurvature"
   return( result )
}
