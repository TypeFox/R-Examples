### 1R
# classification and regression
# continous and discrete data
oneR <- function(formula, data) {
	
	new_data = get.data.frame.from.formula(formula, data)
	new_data = discretize.all(formula,new_data)
	
	class_data = new_data[[1]]
	new_data = new_data[-1] #new_data without class attr
	
	results = sapply(new_data, function(vec) {
			vec = factor(vec)
			errors = sapply(levels(vec), function(val) {
					cvaluestab = as.vector(table(class_data[ which(vec == val) ]))
					return(sum(cvaluestab[ -which.max(cvaluestab) ]) / sum(cvaluestab))
				})
			return(sum(errors))
		})
	
	results = max(results) + min(results) - results
	
	attr_names = dimnames(new_data)[[2]]
	return(data.frame(attr_importance = results, row.names = attr_names))
}
