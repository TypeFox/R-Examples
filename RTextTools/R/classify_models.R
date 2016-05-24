classify_models <- function(container, models, ...) {
# helper method to make it easier to classify with models by algorithm name(s)
# output is a cbinded matrix of model predictions
# hopefully, this method can disappear after refactoring train_model
	result = NULL
	for (name in names(models)) {
		model = models[[name]]
		pred = classify_model(container, model, ...)

		if (is.null(result)) result=pred
		else result = cbind(result, pred)
	}
	
	return(result)
}