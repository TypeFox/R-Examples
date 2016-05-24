itemData.character <-
function(thresholds, p.type = NULL, equation = NULL, ...) {
	model <- CQmodel(show = thresholds, p.type = p.type, equation = equation)
	return(itemData(model,...))
}
