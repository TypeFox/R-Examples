aa.frequency.matrix.function <-
function(aa.count, seqlength){
	position.count = colSums(aa.count)
	
	#Creates a frequency matrix by dividing the number of each amino acid in each position by the total number of nucleotides at each position
	frequency.matrix <- t(aa.count)/position.count
	frequency.matrix <- t(frequency.matrix)
	rownames(frequency.matrix) <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
	return(frequency.matrix)
}
