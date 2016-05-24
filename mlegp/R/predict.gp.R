`predict.gp` <-
function(object, newData = object$X, se.fit = FALSE, ...) {
	if (is.null(dim(newData))) {
		if (!se.fit) {
			return (predictNewZ(object, newData))
		}
		return (list(fit = predictNewZ(object,newData), se.fit = calcPredictionError(object, newData))) 
	}

	ans = matrix(0,dim(newData)[1])
	if (se.fit) se = matrix(0, dim(newData)[1])
	for (i in 1:dim(newData)[1]) {
			ans[i] = predictNewZ(object, newData[i,])
			if (se.fit) se[i] = calcPredictionError(object, newData[i,])
	}
	if (se.fit) return (list(fit = ans, se.fit = se))
	return (ans)
}

