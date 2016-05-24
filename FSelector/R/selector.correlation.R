### correlation.body
# regression
# continous data
correlation.body <- function(formula, data, type = c("pearson", "spearman")) {
	type = match.arg(type)
	
	new_data = get.data.frame.from.formula(formula, data)
	
	lapply(new_data, function(vec) {
		if(is.factor(vec))
			stop("All data must be continous.")
	})
	
	class_data = new_data[[1]]
	new_data = new_data[-1] #new_data without class attr
	
	class_data_complete = complete.cases(class_data)
	results = abs(sapply(new_data, function(attr_data) {
			complete = complete.cases(attr_data) & class_data_complete
			if(!any(complete))
				return(NA)
			vec1 = class_data[complete]
			vec2 = attr_data[complete]
			if(sd(vec1) == 0 || sd(vec2) == 0)
				return(NA)
			return(cor(vec1, vec2, method = type))
		}))
	
	attr_names = dimnames(new_data)[[2]]
	return(data.frame(attr_importance = results, row.names = attr_names))
}

linear.correlation <- function(formula, data) {
	return(correlation.body(formula, data, "pearson"))
}

rank.correlation <- function(formula, data) {
	return(correlation.body(formula, data, "spearman"))
}
