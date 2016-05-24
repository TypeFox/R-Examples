as.simple.formula <- function(attributes, class) {
	return(as.formula(paste(class, paste(attributes, sep = "", collapse = " + "), sep = " ~ ")))
}

get.data.frame.from.formula <- function(formula, data) {
	d = model.frame(formula, data, na.action = NULL)
	for(i in 1:dim(d)[2]) {
		if(is.factor(d[[i]]) || is.logical(d[[i]]) || is.character(d[[i]]))
			d[[i]] = factor(d[[i]])
	}
	return(d)
}

entropyHelper <- function(x) {
    return(entropy(table(x, useNA="always")))
}
