### Consistency
# classification and regression
# continous and discrete data
consistency <- function(formula, data) {
	# uses parent.env
	evaluator <- function(attrs) {
		attr_data = new_data[, attrs, drop=FALSE]

		#deprecated
		#attr_data = as.data.frame(lapply(attr_data, function(vec) { addNA(vec, ifany=TRUE) }))
		#result = aggregate(new_data[[1]], as.list(attr_data), function(classes) {

		hashvec = as.factor(apply(attr_data, 1, digest))
		result = aggregate(new_data[[1]], list(hash=hashvec), function(classes) {
				return(max(as.vector(table(classes))))
			})
		result = sum(result[[dim(result)[2]]]) / dim(attr_data)[1]
		return(result)
	}
	new_data = get.data.frame.from.formula(formula, data)
	new_data = discretize.all(formula, new_data)
	column_names = names(new_data)
	
	result = best.first.search(column_names[-1], evaluator)
	return(result)
}
