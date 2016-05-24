## print a few results of the whole system
print.systemfit <- function( x,
      digits = max( 3, getOption("digits") - 1 ),... ) {

   cat("\n")
   cat("systemfit results \n")
   cat("method: ")
   if(!is.null(x$iter)) if(x$iter>1) cat("iterated ")
   cat( paste( x$method, "\n\n"))
   if(!is.null(x$iter)) {
      if(x$iter>1) {
         if(x$iter<x$control$maxiter) {
            cat( paste( "convergence achieved after",x$iter,"iterations\n\n" ) )
         } else {
            cat( paste( "warning: convergence not achieved after", x$iter,
                        "iterations\n\n" ) )
         }
      }
   }
   cat( "Coefficients:\n" )
   print( x$coefficients, digits = digits )
   invisible( x )
}


## print a few results for a single equation
print.systemfit.equation <- function( x,
      digits = max( 3, getOption("digits") - 1 ), ... ) {

   cat("\n")
   cat( x$method, " estimates for '", x$eqnLabel,
            "' (equation ", x$eqnNo, ")\n", sep = "" )

   cat("Model Formula: ")
   print( formula( x$terms ) )
   if(!is.null(x$inst)) {
      cat("Instruments: ")
      print(x$inst)
   }

   cat("\nCoefficients:")
   print( x$coefficients, digits = digits )
   invisible( x )
}
