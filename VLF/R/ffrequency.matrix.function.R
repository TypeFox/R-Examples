ffrequency.matrix.function <-
function(count.matrix, seqlength){
	#Counts the number of nucleotides in each position of the sequence
	position.count = colSums(count.matrix)
	
	#Creates a frequency matrix by dividing the number of A's, G's, C's, and T's in 	each position by the total number of nucleotides at each position
	frequency.matrix <- t(count.matrix)/position.count
	frequency.matrix <- t(frequency.matrix)
	rownames(frequency.matrix) <- c("A", "C", "G", "T")
	return(frequency.matrix)
}
