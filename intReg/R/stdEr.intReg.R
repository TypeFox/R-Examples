stdEr.intReg <- function( x, boundaries=FALSE, ... ) {
   ## stdEr method.  By default, ignore the fixed boundaries.
   stde <- NextMethod("stdEr", x, ...)
   if(!boundaries) {
      i <- rep(TRUE, nParam(x))
      i[x$param$index$boundary] <- FALSE
      stde <- stde[i]
   }
   class( stde ) <- c( "stdEr.intReg", class(stde) )
   return( stde )
}
