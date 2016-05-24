## gradient function:
## sum over possible individual gradients
logLikGrad <- function(theta, fnOrig, gradOrig, hessOrig,
                       start = NULL, fixed = NULL, sumObs = TRUE,
                       gradAttr = NULL,
                       ...) {

   # Argument "hessOrig" is just for compatibility with logLikHess()
   # argument "gradAttr" should be
   #    - FALSE if the gradient is not provided as attribute of the log-lik value
   #    - TRUE  if the gradient is provided as attribute of the log-lik value
   #    - NULL  if this is not known
   theta <- addFixedPar( theta = theta, start = start, fixed = fixed, ...)

   if(!is.null(gradOrig)) {
      g <- gradOrig(theta, ...)
   } else if( isTRUE( gradAttr ) || is.null( gradAttr ) ) {
      if( exists( "lastFuncGrad" ) && exists( "lastFuncParam" ) ) {
         if( identical( theta, lastFuncParam ) ) {
            g <- lastFuncGrad
         } else {
            g <- "different parameters"
         }
      } else {
         g <- "'lastFuncGrad' or 'lastFuncParam' does not exist"
      }
      if( is.character( g ) ) { # do not call fnOrig() if 'lastFuncGrad' is NULL
         g <- attr( fnOrig( theta, ... ), "gradient" )
      }
   } else {
      g <- NULL
   }
   if( is.null( g ) ) {
      g <- numericGradient(logLikFunc, theta, fnOrig = fnOrig,
         sumObs = sumObs, ...)
   }
   if( sumObs ) {
      ## We were requested a single (summed) gradient.  Return a vector
      g <- sumGradients( g, length( theta ) )
      names( g ) <- names( theta )
      if( !is.null( fixed ) ) {
         g <- g[ !fixed ]
      }
   }
   else {
      ## we were requested individual gradients (if possible).  Ensure g is a matrix
      if(observationGradient(g, length(theta))) {
         ## it was indeed by observations
         g <- as.matrix(g)
         colnames( g ) <- names( theta )
         if( !is.null( fixed ) ) {
            g <- g[ , !fixed ]
         }
      }
      else {
         ## it wasn't
         g <- drop(g)
         names(g) <- names(theta)
         if( !is.null( fixed ) ) {
            g <- g[ !fixed ]
         }
      }         
   }
   return( g )
}
