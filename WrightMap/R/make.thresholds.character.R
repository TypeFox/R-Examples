make.thresholds.character <-
function(item.params,design.matrix="normal",...) {
	#print("character")
	return(make.thresholds(CQmodel(show=item.params),...))
}
