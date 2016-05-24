# draw regression line from model to extremes of fit (J. Fox)

# last modified 2 October 2009 by J. Fox

 
regLine <- function(mod, col=palette()[2], lwd=2, lty=1, ...){
	if(!is.null(class(mod$na.action)) && class(mod$na.action) == "exclude") 
		class(mod$na.action) <-"omit"
	coef <- coefficients(mod)
	if (length(coef) != 2) stop("requires simple linear regression")
	x <- model.matrix(mod)[,2]
	y <- fitted.values(mod)
	min <- which.min(x)
	max <- which.max(x)
	lines(c(x[min], x[max]), c(y[min], y[max]), col=col, lty=lty, lwd=lwd, ...)
}
