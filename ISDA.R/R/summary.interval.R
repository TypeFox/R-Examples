summary.interval <-
function(object,...){
	
	mean = meanInterval(object)
	variance = varianceInterval(object)
	standard_deviance = sdInterval(object)
	median = percentileInterval(object,0.5,"T")
	
	return = data.frame(mean,variance,standard_deviance,median)
	
	## return
	class(object) <- "summary.interval"
	return
}

