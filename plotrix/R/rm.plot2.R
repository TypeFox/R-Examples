
raw.means.plot2 <- function(data, col.id, col.offset, col.x, col.value, fun.aggregate = "mean", ... ) {
	
	if(!is.data.frame(data)) stop("data must be a data.frame")
	
	columns <- c(col.id, col.offset, col.x, col.value)
	
	if (any(!(columns %in% colnames(data)))) stop("column not matching the data")
	
	formula.agg <- as.formula(paste(col.value, "~", col.id, "+", col.offset, "+", col.x))
	
	d.new <- aggregate(formula.agg, data = data, FUN = fun.aggregate)
	raw.means.plot(d.new, col.offset = col.offset, col.x = col.x, col.value = col.value, ...)
}
