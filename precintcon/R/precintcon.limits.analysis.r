#' @export 
precintcon.limits.analysis <- function(...) {
	
	l <- list(...)
	
	if (length(l) > 0) {
		
		result <- lapply(l, FUN = function(x) {
			
			mx <- precintcon.max(x)
			mn <- precintcon.min(x)
			
			return(data.frame(max = mx[1], max.date = mx[2], min = mn[1], min.data = mn[2]))
		})

		dataset <- as.list(match.call()[1:length(l)+1])
		dataset <- do.call(rbind.data.frame, lapply(dataset, FUN = as.character))
		colnames(dataset) <- "dataset"

		
		data <- cbind(
				dataset,
				do.call(rbind.data.frame, result))
		
		return(data)
		
	} else {
		stop("Empty dataset in function precintcon.limits.analysis.")
	}
}