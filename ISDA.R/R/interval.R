interval <-
function(min,max) {
	rval= list(
			minValue = min,
			maxValue = max
	)
	class(rval) = "interval"
	return(rval)
}

