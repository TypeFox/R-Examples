## Return #of observations for models
nObs.maxLik <- function(x, ...) {
   if( is.null( x$gradientObs ) ) {
      stop( "cannot return the number of observations:",
         " please re-run 'maxLik' and",
         " provide a gradient function using argument 'grad' or",
         " (if no gradient function is specified) a log-likelihood function",
         " using argument 'logLik'",
         " that return the gradients or log-likelihood values, respectively,",
         " at each observation" )
   } else if( is.matrix( x$gradientObs ) ) {
      return( nrow( x$gradientObs ) )
   } else {
      stop( "internal error: component 'gradientObs' is not a matrix.",
         " Please contact the developers." )
   }
}
