MODE.second.freq <-
function(freq,seqlength){
	second.large <- c()
	for(i in 1:seqlength){
		freq[which.max(freq[,i]),i] <- NA
		second.large[i] <- max(freq[,i], na.rm = TRUE)
	}
	return(second.large)
}
