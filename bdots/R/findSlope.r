find.slope <- function(time, fixations) {
	quant75 <- (max(fixations) - min(fixations)) * .75 + min(fixations)
	quant25 <- (max(fixations) - min(fixations)) * .25 + min(fixations)
	
	quant75.time <- time[which.min(abs(quant75 - fixations))]
	quant25.time <- time[which.min(abs(quant25 - fixations))]
	
	(quant75 - quant25) / (quant75.time - quant25.time)
}