# The 'stcenter' function is used to center and scale a data frame such that
# its mean is 0 and sd is 1.
# The only difference with the R function 'scale' is that it doesn't center and
# scale column by column, but globally, since all the observations come from
# a single process in the case of space time series.

stcenter <- function(data, center=TRUE, scale=TRUE) {
	
	if (center) {
		# Center data
		attr(data, "mean") <- sum(data) / (nrow(data) * ncol(data))
		data <- data - attr(data, "mean")
	}

	if (scale) {
		# Scale data
		attr(data, "sd") <- sqrt(sum(data^2)/(nrow(data)*ncol(data)-1))
		data <- data / attr(data, "sd")
	}

		return(data)	
		
}

