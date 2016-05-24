VLF.aminoAcids <-
function(convert.matrix, seq.matrix, seqlength){
	aa_matrix <- matrix(NA, nrow = nrow(convert.matrix), ncol = seqlength + 2)
	aa_matrix[,1:2] <- seq.matrix[,1:2]
	for(i in 1:seqlength){
		aa_matrix[which(!is.na(convert.matrix[,i+2])),i+2] <- seq.matrix[which(!is.na(convert.matrix[,i+2])),i+2]
	}
	return(aa_matrix)
}
