itemData.default <-
function(thresholds, item.type = "deltas",...) {
	if(item.type == "thresholds")
		thresholds <- make.thresholds(thresholds,...)
	return(thresholds)
}
