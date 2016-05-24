"plot.preplot.gam" <-
function(x, y = NULL, residuals = NULL, rugplot = TRUE, se = FALSE, scale = 0, fit = TRUE,
	...)
{
	listof <- inherits(x[[1]], "preplot.gam")
	if(listof) {
		TT <- names(x)
		scales <- rep(0, length(TT))
		names(scales) <- TT
		for(i in TT)
			scales[i] <- plot.preplot.gam(x[[i]], y = NULL, 
				residuals, rugplot, se, scale, fit, ...)
		#			scales[i] <- UseMethod("plot",x[[i]])
		invisible(scales)
	}
	else {
		dummy <- function(residuals = NULL, rugplot = TRUE, se = FALSE, scale
			 = 0, fit = TRUE, ...)
		c(list(residuals = residuals, rugplot = rugplot, se = se, scale
			 = scale, fit = fit), list(...))
		d <- dummy(residuals, rugplot, se, scale, fit, ...)
		uniq.comps <- unique(c(names(x), names(d)))
		Call <- c(as.name("gplot"), c(d, x)[uniq.comps])
		mode(Call) <- "call"
		invisible(eval(Call))
	}
}
