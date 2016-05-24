specimen.frequencies <-
function(freq, seq.matrix, no.spec, spec.names,seqlength){
	sequence.freq.matrix <- matrix(NA, nrow = no.spec, ncol = seqlength)
	letter_vector<-c("A","C","G","T")
	for(i in 1:seqlength){
		sequence.freq.matrix[which(seq.matrix[,i+2] == letter_vector[1]),i] <- freq[1,i]
		sequence.freq.matrix[which(seq.matrix[,i+2] == letter_vector[2]),i] <- freq[2,i]
		sequence.freq.matrix[which(seq.matrix[,i+2] == letter_vector[3]),i] <- freq[3,i]
		sequence.freq.matrix[which(seq.matrix[,i+2] == letter_vector[4]),i] <- freq[4,i]
	}
	rownames(sequence.freq.matrix) <- spec.names
	return(sequence.freq.matrix)
}
