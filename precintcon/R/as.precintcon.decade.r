#' @export
as.precintcon.decade <- function(object) {

  object <- as.precintcon.annual(object)
	
	result <- data.frame()
	
	start <- head(which(head(object$year, n = 10) %% 10 == 0), n = 1)
	end   <- nrow(object) - tail(which(rev(tail(object$year, n = 10)) %% 10 == 9), n = 1) + 1
	
	print(seq(start, end, by = 10))
	for (i in seq(start, end, by = 10))
		result <- rbind(result, data.frame(year=object[i,1], precipitation=sum(object[i:(i+9),2])))
	
	class(result) <- c("data.frame", "precintcon.decade")
	
	return(result)
}
