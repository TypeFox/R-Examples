### methods for extracting standard errors from the models

stdEr <- function(x, ...)
    ## Extract standard deviations from models (as coefficients)
    UseMethod("stdEr")

stdEr.default <- function(x, ...) {
   if( !isS4( x ) ) {
      if( !is.null( x$std ) ) {
         return(x$std)
      }
   }
   if(!is.null(vc <- vcov(x))) {
      s <- sqrt(diag(vc))
      names(s) <- names(coef(x))
      return(s)
   }
   return(NULL)
                           # if neither std nor vcov is defined, we return NULL...
}

stdEr.lm <- function(x, ...)
    sqrt(diag(vcov(x)))

