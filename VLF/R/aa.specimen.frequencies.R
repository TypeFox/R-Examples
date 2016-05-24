aa.specimen.frequencies <-
function(freq, seq.matrix, spec.names, seqlength){	
	no.spec <- nrow(seq.matrix)
	letter_vector <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
	aasequence.freq.matrix <- matrix(NA, nrow = no.spec, ncol = seqlength)
		for(i in 1:seqlength){
		aasequence.freq.matrix[which(seq.matrix[,i+2] == letter_vector[1]),i] <- freq[1,i]
		aasequence.freq.matrix[which(seq.matrix[,i+2] == letter_vector[2]),i] <- freq[2,i]
		aasequence.freq.matrix[which(seq.matrix[,i+2] == letter_vector[3]),i] <- freq[3,i]
		aasequence.freq.matrix[which(seq.matrix[,i+2] == letter_vector[4]),i] <- freq[4,i]
		aasequence.freq.matrix[which(seq.matrix[,i+2] == letter_vector[5]),i] <- freq[5,i]
		aasequence.freq.matrix[which(seq.matrix[,i+2] == letter_vector[6]),i] <- freq[6,i]
		aasequence.freq.matrix[which(seq.matrix[,i+2] == letter_vector[7]),i] <- freq[7,i]
		aasequence.freq.matrix[which(seq.matrix[,i+2] == letter_vector[8]),i] <- freq[8,i]
		aasequence.freq.matrix[which(seq.matrix[,i+2] == letter_vector[9]),i] <- freq[9,i]
		aasequence.freq.matrix[which(seq.matrix[,i+2] == letter_vector[10]),i] <- freq[10,i]
		aasequence.freq.matrix[which(seq.matrix[,i+2] == letter_vector[11]),i] <- freq[11,i]
		aasequence.freq.matrix[which(seq.matrix[,i+2] == letter_vector[12]),i] <- freq[12,i]
		aasequence.freq.matrix[which(seq.matrix[,i+2] == letter_vector[13]),i] <- freq[13,i]
		aasequence.freq.matrix[which(seq.matrix[,i+2] == letter_vector[14]),i] <- freq[14,i]
		aasequence.freq.matrix[which(seq.matrix[,i+2] == letter_vector[15]),i] <- freq[15,i]
		aasequence.freq.matrix[which(seq.matrix[,i+2] == letter_vector[16]),i] <- freq[16,i]
		aasequence.freq.matrix[which(seq.matrix[,i+2] == letter_vector[17]),i] <- freq[17,i]
		aasequence.freq.matrix[which(seq.matrix[,i+2] == letter_vector[18]),i] <- freq[18,i]
		aasequence.freq.matrix[which(seq.matrix[,i+2] == letter_vector[19]),i] <- freq[19,i]
		aasequence.freq.matrix[which(seq.matrix[,i+2] == letter_vector[20]),i] <- freq[20,i]		
	}
	rownames(aasequence.freq.matrix) <- spec.names
	return(aasequence.freq.matrix)
}
