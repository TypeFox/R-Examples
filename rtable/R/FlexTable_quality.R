#' @importFrom stats AIC
#' @importFrom stats BIC
#' @importFrom stats logLik
#' @importFrom stats fitted
#' @importFrom stats lm
#' @importFrom stats pf
#' @title get quality measures from a statistical models
#'
#' @description
#' Extract quality measures from a statistical models and return them as a
#' \code{\link{FlexTable}}.
#'
#' @param object a model (supported model are lm, aov, lme, glm, gls, clm, lmerMod, glmerMod, multinom, clmm)
#' @param ... further arguments, not used.
#' @return a \code{FlexTable} object
#' @export
FlexTable_quality = function (object, ...){

	if( inherits(object, c( "lme", "clmm", "glm", "gls", "clm", "lmerMod", "glmerMod", "multinom") ) ){
		quality_data = data.frame( label = c("Observations", "Akaike Inf. Crit.", "Bayesian Inf. Crit.", "Log Likelihood"),
				value =
						c( as.character( length( fitted(object) ) ),
						format(AIC(object)),
						format( BIC(object) ) ,
						format( logLik(object) ) )
		)
	} else if( inherits(object$fit, c( "lm", "aov") ) ){
		quality_data = get_mod_lin_quality( object )
	} else stop("unimplemented model")

	ft = FlexTable( quality_data, header.columns = F )
	ft[,1] = parRight(padding=2)
	ft[,1] = textItalic()
	ft[,2] = parLeft(padding=2)

	ft = setFlexTableBorders(ft, footer = T,
		inner.vertical = borderNone(), inner.horizontal = borderNone(),
		outer.horizontal = borderNone(), outer.vertical = borderNone())

	ft
}

get_mod_lin_quality = function( x, digits = max(3L, getOption("digits") - 3L) ){
  if( inherits(x, "aov" ) )
    x = lm( x )
  if( !inherits(x, "lm" ) )
    stop("not a linear model")

  out = matrix( c("Observations", as.character( length( x$"fitted.values" ) ) ), ncol= 2 )

  x = summary(x)

  if (!is.null(x$fstatistic)) {
    fstat = paste( formatC(x$fstatistic[1L], digits = digits), "on",
      x$fstatistic[2L], "and", x$fstatistic[3L], "DF, p-value:",
  	    format.pval(pf(x$fstatistic[1L], x$fstatistic[2L],
          x$fstatistic[3L], lower.tail = FALSE), digits = digits ))
    qual_crit = list(
        c("Multiple R-squared", formatC(x$r.squared, digits = digits) ),
        c("Adjusted R-squared", formatC(x$adj.r.squared, digits = digits) ),
        c("F-statistic", fstat )
    )
    qual_crit = do.call( rbind, qual_crit)
    out = rbind( out, qual_crit )
  }
  out = as.data.frame( out, stringsAsFactors = F )
  names( out ) = c("label", "value")
  out
}
