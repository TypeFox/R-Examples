"gplot.default" <-
function(x, y, se.y = NULL, xlab = "", ylab = "", residuals = NULL, rugplot = FALSE,
	scale = 0, se = FALSE, fit = TRUE, ...)
switch(data.class(x)[1],
       AsIs = { class(x)<-NULL
                gplot.default(x , y = y, se.y = se.y, xlab = xlab,
		ylab = ylab, residuals = residuals, rugplot = rugplot, scale = 
		scale, se = se, fit = fit, ...)
              },
	logical = gplot.factor(x = factor(x), y = y, se.y = se.y, xlab = xlab,
		ylab = ylab, residuals = residuals, rugplot = rugplot, scale = 
		scale, se = se, fit = fit, ...),
	list = gplot.list(x = x, y = y, se.y = se.y, xlab = xlab, ylab = ylab,
		residuals = residuals, rugplot = rugplot, scale = scale, se = 
		se, fit = fit, ...),
	if(is.numeric(x)) gplot.numeric(x = as.vector(x), y = y, se.y = se.y,
			xlab = xlab, ylab = ylab, residuals = residuals, 
			rugplot = rugplot, scale = scale, se = se, fit = fit,
			...) else warning(paste("The \"x\" component of \"",
			ylab, "\" has class \"", paste(class(x), collapse = 
			"\", \""), "\"; no gplot() methods available", sep = ""
			)))
 
