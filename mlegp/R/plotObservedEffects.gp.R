`plotObservedEffects.gp` <-
function(x, params = 1:x$numDim, ...) {
	createWindow(length(params))
	params = toParamIndexes(params, x$params)
	for (i in params) {
		plot(x$X[,i], x$Z, xlab = x$params[i])
	}
}

