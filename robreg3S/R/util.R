
print.robreg3S <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("Coefficients:\n")
  print(format(x$Summary.Table, digits = digits), print.gap = 2, quote = FALSE)
  cat("\n")
  invisible(x)
}

coef.robreg3S <- function(object, ...) object$coef

summary.robreg3S <- function(object, ...) object$Summary.Table

confint.robreg3S <- function(object, ...) {
	if( all(is.na(object$Summary.Table[,2] )) ) return( NA )
	cbind( object$Summary.Table[,1] - qnorm(0.975) * object$Summary.Table[,2],
		object$Summary.Table[,1] + qnorm(0.975) * object$Summary.Table[,2] )
	
}