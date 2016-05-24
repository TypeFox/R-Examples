preProcess.removeBias = function(object,...) {
    params = object$params	
	if (!is.na(params$removeBias[[1]])) 
		object = biasCorr(object,...)
	object
}
