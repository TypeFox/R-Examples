coef.sma <- function(object, ...){

	x <- object
	if(length(x$coef) == 1){
		res <- x$coef[[1]][,1]
		names(res) <- c("elevation","slope")
	} else {
		res <- lapply(x$coef, "[", 1)
		res <- as.data.frame(do.call("rbind",lapply(res, t)))
		rownames(res) <- names(x$coef)
	}
	return(res)
	
}