Declutter <- function(tmp, delta) {
	# function that takes data from early experiments where double points are recorded
	if (missing(delta)) 
		delta <- 2/60
	# button presses within 2s of each other
	
	if (class(tmp) == "dark") {
		if (is.null(tmp$data)) 
			tmp$data = "unknown"

		# idx <- c(diff(tmp$time),1) > delta    
		# this removes the first button press

		idx <- c(1, diff(tmp$time)) > delta
		# this removes the second press
		
		x <- tmp$time[idx]
		y <- tmp$thrs[idx]

	} else {
		idx <- c(0, diff(tmp[, 1])) > delta
		x <- tmp[idx, 1]
		y <- tmp[idx, 2]

	}
		tmp$time = x
		tmp$thrs = y
		obj <- tmp
		class(obj) = "dark"
		return(obj)
}