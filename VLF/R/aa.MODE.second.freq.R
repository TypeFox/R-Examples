aa.MODE.second.freq <-
function(freq.matrix, seqlength){
	second.large <- c()
	for(i in 1:seqlength){
		freq.matrix[which.max(freq.matrix[,i]),i] <- NA
		second.large[i] <- max(freq.matrix[,i], na.rm = TRUE)
	}
	return(second.large)
}
