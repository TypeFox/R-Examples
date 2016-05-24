MODE <-
function(freq, seqlength){
	sequence.ind <- c("A", "C", "G", "T")
	biggest.sequence.num <- c()
	biggest.sequence.nuc <- c()
	for(i in 1:seqlength){
		biggest.sequence.num[i] <- which.max(freq[,i])
		biggest.sequence.nuc[i] <- sequence.ind[biggest.sequence.num[i]]
	}
	return(biggest.sequence.nuc)
}
