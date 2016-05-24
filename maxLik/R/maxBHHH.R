maxBHHH <- function(fn, grad=NULL, hess=NULL,
                    start,
                    finalHessian="BHHH",
                    ...) {
   ## hess:   Hessian, not used, for compatibility with the other methods

   ## check if arguments of user-provided functions have reserved names
   argNames <- c( "fn", "grad", "hess", "start", "print.level", "iterlim" )
   checkFuncArgs( fn, argNames, "fn", "maxBHHH" )
   if( !is.null( grad ) ) {
      checkFuncArgs( grad, argNames, "grad", "maxBHHH" )
   }
   if( !is.null( hess ) ) {
      checkFuncArgs( hess, argNames, "hess", "maxBHHH" )
   }

   ## using the Newton-Raphson algorithm with BHHH method for Hessian
   a <- maxNR( fn=fn, grad = grad, hess = hess, start=start,
              finalHessian = finalHessian,
              bhhhHessian = TRUE,
               ...)

   a$type = "BHHH maximisation"
   invisible(a)
}
