aa.MODE <-
function(freq.matrix, seqlength){
	sequence.ind <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
	biggest.sequence.num <- c()
	biggest.sequence.aa <- c()
	for(i in 1:seqlength){
		biggest.sequence.num[i] <- which.max(freq.matrix[,i])
		biggest.sequence.aa[i] <- sequence.ind[biggest.sequence.num[i]]
	}
	return(biggest.sequence.aa)
}
