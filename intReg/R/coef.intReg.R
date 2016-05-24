coef.intReg <- function( object, boundaries=FALSE, ... ) {
   ## coef method.  By default, ignore the fixed boundaries.
   coefValues <- NextMethod("coef", object)
   if(is.null(coefValues))
                           # may be NULL in case of certain errors
                           # (or somebody hacking the variables)
                           # (should add a test)
       return(coefValues)
   if(!boundaries) {
      i <- rep(TRUE, length(coefValues))
      i[object$param$index$boundary] <- FALSE
      coefValues <- coefValues[i]
   }
   class( coefValues ) <- c( "coef.intReg", class(coefValues) )
   return( coefValues )
}
