## calculate confidence intervals of the coefficients
confint.systemfit <- function( object, parm = NULL, level = 0.95,
      useDfSys = NULL, ... ) {

   if( is.null( useDfSys ) ) {
      useDfSys <- length( coef( object ) ) != object$rank
         # TRUE if there are restrictions imposed
   }

   probLower <- ( 1 - level ) / 2
   probBoth <- c( probLower, 1 - probLower )
   pct <- paste( round( 100 * probBoth, 1 ), "%" )
   ci <- matrix( NA, length( object$coefficients ), 2,
            dimnames = list( names( object$coefficients ), pct ) )
   j <- 1
   for( i in 1:length( object$eq ) ) {
      object$eq[[i]]$df.residual.sys <- object$df.residual
      ci[ j:(j+length( coef( object$eq[[ i ]] ) )-1), ] <- confint( object$eq[[ i ]],
         useDfSys = useDfSys )
      j <- j + length( coef( object$eq[[ i ]] ) )
   }
   class( ci ) <- "confint.systemfit"
   ci
}

## calculate confidence intervals of the coefficients of a single equation
confint.systemfit.equation <- function( object, parm = NULL, level = 0.95,
   useDfSys = NULL, ... ) {

   if( is.null( useDfSys ) ) {
      useDfSys <- object$nCoef.sys != object$rank.sys
         # TRUE if there are restrictions imposed
   }

   probLower <- ( 1 - level ) / 2
   probBoth <- c( probLower, 1 - probLower )
   pct <- paste( round( 100 * probBoth, 1 ), "%" )
   ci <- matrix( NA, length( object$coefficients ), 2,
            dimnames = list( names( object$coefficients ), pct ) )
   if( useDfSys ) {
      fac <- qt( probBoth, object$df.residual.sys )
   } else {
      fac <- qt( probBoth, object$df.residual )
   }
   coef <- summary( object )$coefficients
   ci[] <- coef[ , 1 ] + coef[ , 2 ] %o% fac
   class( ci ) <- "confint.systemfit"
   ci
}

## print the confidence intervals of the coefficients
print.confint.systemfit <- function( x, digits = 3, ... ) {
   print( unclass( round( x, digits = digits, ...) ) )
   invisible(x)
}

