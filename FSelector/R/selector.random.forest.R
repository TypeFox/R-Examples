### RANDOM FOREST
# classification and regression
# continous and discrete data
# NA deleted
random.forest.importance <- function(formula, data, importance.type = 1) {
	new_data = get.data.frame.from.formula(formula, data)

	# get rid of NAs
	no_na = rep(TRUE, dim(new_data)[1])
	for(i in 1:dim(new_data)[2]) {
		no_na = no_na & complete.cases(new_data[, i])
	}
	new_data = new_data[no_na, , drop=FALSE]
	
	forest = randomForest(formula, new_data,
		ntree = 1000, keep.forest = FALSE, importance = TRUE)
		
	res = as.data.frame(importance(forest, type = importance.type))
	colnames(res)[1] = "attr_importance"
	return(res)
}
