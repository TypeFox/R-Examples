# quant calculates the posterior quantile of the true value of theta
quant <- function (draws) {
	n<-length(draws)
	rank.theta <- c(1:n)[order(draws)==1] - 1
	quantile.theta <- (rank.theta +.5) / n
	return(quantile.theta) }

