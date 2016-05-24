vcov.intReg <- function( object, boundaries=FALSE, ... ) {
   ## vcov method.  By default, ignore the fixed boundaries.
   vc <- NextMethod("vcov", object, ...)
   if(!boundaries) {
      i <- rep(TRUE, nParam(object))
      i[object$param$index$boundary] <- FALSE
      vc <- vc[i,i]
   }
   class( vc) <- c( "vcov.intReg", class(vc) )
   return( vc )
}
