summary.bezierCurveFit <- function(object, ...){
	r <- ''

	r <- c(r, '\nbezierCurveFit Summary\n')
	
	for(i in 1:length(object$p)){
		r <- c(r, '\tDimension ', i, ' fit with ', length(object$p[[i]]), ' parameters and RSE of ', format(object$rse[[i]]), '.')
		r <- c(r, ' ', object$fit.stopped.by[i], '.')
		r <- c(r, '\n')
	}

	r <- c(r, '\n')
	class(r) <- "summary.bezierArcLength"
	r
}

print.summary.bezierCurveFit <- function(x, ...) cat(x, sep='')