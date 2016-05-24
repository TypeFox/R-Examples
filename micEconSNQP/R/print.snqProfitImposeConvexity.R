print.snqProfitImposeConvexity <- function( x, digits=6,... ) {

   save.digits <- unlist( options( digits = digits ) )
   on.exit( options( digits = save.digits ) )

   table <- NULL
   labels <- NULL

   cat( "\nEstimation results of the SNQ profit function with convexity" )
   cat( "imposed:\n")
   cat( "Functional form: " )
   if( x$form == 0 ) {
      cat( "default (" )
   } else {
      cat( "extended (" )
   }
   cat( x$form, ")\n", sep = "" )
   cat( "Minimum distance estimation:\n" )
   if( x$mindist$convergence == 0 ) {
      cat( paste( "convergence achieved after", x$mindist$count[ 1 ],
         "iterations\n" ) )
   } else {
      cat( paste( "warning: convergence not achieved after (code ",
         x$mindistest$convergence, ")\n" ) )
   }
   cat( "\nEstimated coefficients:\n" )
   coef <- matrix( x$coef$allCoef, ncol = 1 )
   rownames( coef ) <- names( x$coef$allCoef )
   print( coef )

   cat( "\nR-squared values of the netput equations:\n" )
   print( x$r2 )

   cat( "\nPrice elasticities of the netputs:\n" )
   print( x$ela$ela )

   if( x$convexity ) {
      #cat( "\nThis profit function is convex in netputs." )
   } else {
      cat( "\nwarning: this profit function is not convex in netput prices.\n" )
   }
   cat( "\n" )

   invisible( x )
}
