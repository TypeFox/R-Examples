#' @importFrom stats coef
#' @title get coefficients from a statistical models
#'
#' @description
#' Extract coefficients from a statistical models and return them as a
#' \code{\link{FlexTable}}.
#'
#' @param object a model (supported model are lm, aov, lme, glm, gls, clm, lmerMod, glmerMod, multinom, clmm)
#' @param ... further arguments, not used.
#' @return a \code{FlexTable} object
#' @export
FlexTable_coef = function (object, ...){

	out = NULL
	object      = object
	digits = NULL
	if( inherits( object, c("glmerMod", "lmerMod") ) ){
  	out = coef( summary(object) )
		digits = c(0, 3, 3, 3)
	} else if( inherits(object, c("aov", "lm", "glm", "clm", "clmm" ))) {
		out = coef( summary(object) )
		digits = c(0, 3, 3, 3, 3)
	} else if( inherits(object, "lme") ) {
		out = summary(object)$tTable
		digits = c(0, 3, 3, 0, 3, 3)
	} else if( inherits(object, "gls" ) ) {
		out = summary(object)$tTable
		digits = c(0, 3, 3, 3, 3)
	} else if( inherits(object, "multinom")) {
		x = summary( object )
		if(x$is.binomial) {
			out = cbind(Values = x$coefficients,
					"Std. Err." = x$standard.errors,
					"Value/SE" = x$Wald.ratios)
			digits = c(0, 3, 3 )
		} else {
			out = cbind(Values = x$coefficients,
					"Std. Err." = x$standard.errors )
			digits = c( 0, 3 )
			if(!is.null(x$Wald.ratios)) {
				out = cbind( out, "Value/SE" = x$coefficients/x$standard.errors )
				digits = c(digits, 3 )
			}
		}
	}
	else stop("unimplemented object")

	ttable= xtable(out)
	if( !is.null( digits ))
		digits(ttable) = digits
	as.FlexTable( ttable )
}
