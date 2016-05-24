VLF.convert.matrix <-
function(seq.matrix, freq, p,seqlength){
	converted <- matrix(NA, nrow = nrow(freq), ncol = seqlength+2)
	converted[,1:2] = seq.matrix[,1:2]
	for(i in 1:seqlength){
		converted[which(freq[,i] < p), i+2] <- freq[which(freq[,i] < p), i]
	}
	return(converted)
}
