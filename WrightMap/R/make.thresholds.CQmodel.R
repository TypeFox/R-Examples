make.thresholds.CQmodel <-
function(item.params,item.table = NULL,interactions = NULL,step.table = NULL,design.matrix = "normal",throld = .5,alpha = 1,...) {
	deltas <- make.deltas(item.params,item.table,interactions,step.table)
	#message("Creating threshold parameters out of deltas")
	return(make.thresholds(deltas,throld = throld,alpha = alpha,...))
}
