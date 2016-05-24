MODE.freq <-
function(freq,seqlength){
	largest.freq <- c()
	for(i in 1:seqlength){
		largest.freq[i] <- max(freq[,i])
	}
	return(largest.freq)
}
