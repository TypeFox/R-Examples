estfun.maxLik <- function( x, ... ) {

   if( is.null( x$gradientObs ) ) {
      stop( "cannot return the gradients of the log-likelihood function",
         " evaluated at each observation: please re-run 'maxLik' and",
         " provide a gradient function using argument 'grad' or",
         " (if no gradient function is specified) a log-likelihood function",
         " using argument 'logLik'",
         " that return the gradients or log-likelihood values, respectively,",
         " at each observation" )
   }

   return( x$gradientObs )
}
