calculate.sample.prev <-
function(data, status, tag.healthy) {
	sample.prev <- nrow(data[data[,status] != tag.healthy,])/nrow(data)
  	if (sample.prev == 0) {
		stop("There are no diseased subjects in your dataset. Please review data and \n variables. Prevalence must be a value higher than 0 and lower than 1.", call. = FALSE)	 		
  	}
  	if (sample.prev == 1) {
		stop("There are no healthy subjects in your dataset. Please review data and \n variables. Prevalence must be a value higher than 0 and lower than 1.", call. = FALSE)			
  	}
  	res <- sample.prev
}
