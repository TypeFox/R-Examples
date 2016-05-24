coef.aidsEst <- function( object, ... ) {
   result <- list()
   result$alpha0 <- object$coef$alpha0
   result$alpha  <- object$coef$alpha
   result$beta   <- object$coef$beta
   result$gamma  <- object$coef$gamma
   result$delta  <- object$coef$delta
   class( result ) <- "coef.aidsEst"
   return( result )
}
