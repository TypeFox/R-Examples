calculate.empirical.AUC <-
function(data, marker, status, tag.healthy, direction = c("<",">"), conf.level = 0.95) {
	direction <- match.arg(direction)
	
	marker.diseased = data[data[,status] != tag.healthy, marker]
	n.diseased = length (marker.diseased)

	marker.healthy = data[data[,status] == tag.healthy, marker]
	n.healthy = length(marker.healthy)

	d.h <- sapply(marker.healthy, function(x, m.d, direction = c("<",">")) {
		diff <- outer(x, m.d, "-")
		res <- vector(length = 2)
		res[1] <- if(direction == "<") {
			sum(diff < 0)
		} else {
			sum(diff > 0)
		}
		res[2] <- sum(diff == 0)
		res
	}, m.d = marker.diseased, direction = direction)
	
	d.d <- sapply(marker.diseased, function(x, m.h, direction = c("<",">")) {
		diff <- - outer(x, m.h, "-")
		res <- vector(length = 2)
		res[1] <- if(direction == "<") {
			sum(diff < 0)
		} else {
			sum(diff > 0)
		}
		res[2] <- sum(diff == 0)
		res
	}, m.h = marker.healthy, direction = direction)
	
	area <- (sum(d.h[1,]) + 0.5*sum(d.h[2,]))/(n.diseased*n.healthy)
	sum1 <- sum(((d.d[1,]+0.5*d.d[2,])/n.healthy-area)^2)
	first.term <- sum1/(n.diseased*(n.diseased-1))
	sum2 <- sum(((d.h[1,]+0.5*d.h[2,])/n.diseased-area)^2)
	second.term <- sum2/(n.healthy*(n.healthy-1))
	var <- first.term+second.term
	z <- qnorm(1-((1-conf.level)/2))	 	
	# Lower end of (1-conf.level)% confidence interval:		 
	inf <- area-z*sqrt(var)
	# Upper end of (1-conf.level)% confidence interval:	 	
	sup <- area+z*sqrt(var)

	res <- c(area, inf, sup)
	names(res) <- c("AUC", "ll", "ul")
	return(res)
}
