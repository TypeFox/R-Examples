#' @export
precintcon.rai.analysis <- function(object, granularity = "m") {
	
	tmp <- NULL
	
	if (granularity == "m") {
		object <- as.precintcon.monthly(object);
	} else if (granularity == "a") {
		object <- as.precintcon.annual(object);
	} else {
		stop("granularity should be either 'm' or 'a'.")
	}

	object[is.na(object)] <- 0
	
	tmp <- sort(object$precipitation)
	tmp <- tmp[tmp >= 0]
	
	l <- length(tmp)
	
	xb <- mean(tmp[1:10])
	mb <- mean(tmp[(l-9):l])
	k  <- mean(object$precipitation)

	result <- apply(data.frame(object$precipitation), 1, FUN = function(x) {
	
		result <- 
			if (x < k)
				-3 * ((x - k) / (xb - k))
			else
				 3 * ((x - k) / (mb - k))		
	
		return(result)
		
	})
	
	result <- cbind(object[1:(ncol(object)-1)], data.frame(rai = result))

	class(result) <- c("data.frame", "precintcon.rai")
	
	return(result)
}