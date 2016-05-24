`SwapOrder` <-
function(data) {
	d <- data
	if (inherits(d, "mlds.df")) {
		invord <- attr(data, "invord")
		} else
	{
	invord <- (data$S1 > data$S4) }
	d$resp[invord] <- 1 - data$resp[invord]
	d[invord, -1] <- d[invord, c(4, 5, 2, 3)]
	d	
}

