normalize.min.max <- function(data) {
	attr_count = dim(data)[2]
	if(attr_count == 0)
		return(data)
	for(i in 1:attr_count) {
		if(!is.numeric(data[, i]))
			next()
		if(!any(complete.cases(data[, i])))
			next()
		mm = range(data[, i], na.rm = TRUE)
		minimum = mm[1]
		maximum = mm[2]
		if(minimum == maximum)
			data[, i] = data[, i] / minimum
		else
			data[, i] = (data[, i] - minimum) / (maximum - minimum)
	}

	return(data)
}
