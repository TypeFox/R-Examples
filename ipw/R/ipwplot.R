ipwplot <- function(weights, timevar = NULL, binwidth = NULL, logscale = TRUE, xlab = NULL, ylab = NULL, main = "", ref = TRUE, ...){
	if (!is.null(timevar)){
		timevargrp <- round(timevar/binwidth)*binwidth
		if (is.null(xlab)) xlab <- deparse(match.call()$timevar, width.cutoff = 500)
		if (is.null(ylab) & logscale == FALSE) ylab <- deparse(match.call()$weights, width.cutoff = 500)
		if (is.null(ylab) & logscale == TRUE) ylab <- paste("log(", deparse(match.call()$weights, width.cutoff = 500), ")", sep = "")
		if(logscale == TRUE){
			boxplot(log(weights) ~ timevargrp, pch = 20, xlab = xlab, ylab = ylab, main = main, ...)
			if (ref == TRUE) abline(h = log(1), lty = 2)
		}
		if(logscale == FALSE){
			boxplot(weights ~ timevargrp, pch = 20, xlab = xlab, ylab = ylab, main = main, ...)
			if (ref == TRUE) abline(h = 1, lty = 2)
		}
	}
	if (is.null(timevar)){
		if (is.null(xlab) & logscale == FALSE) xlab <- deparse(match.call()$weights, width.cutoff = 500)
		if (is.null(xlab) & logscale == TRUE) xlab <- paste("log(", deparse(match.call()$weights, width.cutoff = 500), ")", sep = "")
		if (is.null(ylab)) ylab <- "density"
		if(logscale == TRUE){
			plot(density(log(weights)), pch = 20, xlab = xlab, ylab = ylab, main = main, ...)
			if (ref == TRUE) abline(v = log(1), lty = 2)
		}
		if(logscale == FALSE){
			plot(density(weights), pch = 20, xlab = xlab, ylab = ylab, main = main, ...)
			if (ref == TRUE) abline(v = 1, lty = 2)
		}
	}
}
