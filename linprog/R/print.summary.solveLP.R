## print the summary results
print.summary.solveLP <- function( x, ... ) {
   object <- x

   cat("\n\nResults of Linear Programming / Linear Optimization\n")

   cat("\nObjective function")
   if( object$maximum ) {
      cat(" (Maximum): ")
   } else {
      cat(" (Minimum): ")
   }
   cat( object$opt, "\n" )

   cat("\nSolution\n")
   object$solution <- as.matrix(object$solution)
   colnames( object$solution ) <- c("opt")
   print( object$solution )
   cat("\n")
   invisible( x )
}
