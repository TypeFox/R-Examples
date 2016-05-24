ts2df <-
function(x, mon1 = 1, addYr = FALSE, omit = FALSE) {

	# validate args
	if (!is(x, 'ts') || is(x, 'mts') || !identical(frequency(x), 12))
		stop("x must be a monthly 'ts' vector")
	if (!mon1 %in% 1:12)
		stop("mon1 must be between 1 and 12")
		
	# convert to data.frame
	x1 <- window(x, start = c(start(x)[1] - 1, mon1), end = c(end(x)[1] +
		1, ifelse(mon1 == 1, 12, mon1 - 1)), extend = TRUE)
	d1 <- as.data.frame(matrix(x1, byrow = TRUE, ncol = 12))	
	colnames(d1) <- if (mon1 == 1) month.abb else month.abb[c(mon1:12,
		1:(mon1 - 1))]
	rownames(d1) <- (start(x1)[1] + addYr):(start(x1)[1] + nrow(d1) - 1 +
		addYr)
	
	# trim leading and trailing NA rows, and optionally, rows with any NAs
	d1 = d1[apply(d1, 1, function(x) !all(is.na(x))),]
	if (omit) d1 = na.omit(d1)
	d1
}
