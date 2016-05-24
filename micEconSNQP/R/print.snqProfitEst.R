print.snqProfitEst <- function( x, ... ) {

   table <- NULL
   labels <- NULL

   cat( "\nEstimation results of the SNQ profit function:\n")
   cat( "Functional form: " )
   if( x$form == 0 ) {
      cat( "default (" )
   } else {
      cat( "extended (" )
   }
   cat( x$form, ")\n", sep = "" )
   cat( "Estimation method: " )
   if( !is.null( x$est$iter ) ) {
      if( x$est$iter > 1 ) cat( "iterated " )
   }
   cat( paste( x$est$method, "\n", sep = "" ) )
   if( !is.null( x$est$iter ) ) {
      if( x$est$iter > 1 ) {
         if( x$est$iter < x$est$maxiter ) {
            cat( paste( "convergence achieved after", x$est$iter,
               "iterations\n" ) )
         } else {
            cat( paste( "warning: convergence not achieved after", x$est$iter,
               "iterations\n" ) )
         }
      }
   }
   cat( "\nEstimated coefficients:\n" )
   print( x$coef$stats, ... )

   cat( "\nR-squared values of the netput equations:\n" )
   print( x$r2, ... )

   cat( "\nPrice elasticities of the netputs:\n" )
   print( x$ela$ela, ... )

   if( x$convexity ) {
      #cat( "\nThis profit function is convex in netputs." )
   } else {
      cat( "\nwarning: this profit function is not convex in netput prices.\n" )
   }
   cat( "\n" )

   invisible( x )
}
