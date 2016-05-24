VLF.nucleotides <-
function(convert.matrix, seq.matrix,seqlength){
	nuc_matrix <- matrix(NA, nrow = nrow(convert.matrix), ncol = seqlength+2)
	nuc_matrix[,1:2] <- seq.matrix[,1:2]
	for(i in 1:seqlength){
		nuc_matrix[which(!is.na(convert.matrix[,i+2])),i+2] <- seq.matrix[which(!is.na(convert.matrix[,i+2])),i+2]
	}
	return(nuc_matrix)
}
