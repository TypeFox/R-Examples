aa.MODE.freq <-
function(freq.matrix, seqlength){
	largest.freq <- c()
	for(i in 1:seqlength){
		largest.freq[i] <- max(freq.matrix[,i])
	}
	return(largest.freq)
}
