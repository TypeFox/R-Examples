count.function <-
function(nucleotides, spec.no, seqlength){
	count <- mat.or.vec(nr = 4, nc = seqlength)
	letter_vector <- c("A","C","G","T")
	for(i in 1:seqlength){
		count[1,i] = length(which(nucleotides[,i+2] == letter_vector[1]))
		count[2,i] = length(which(nucleotides[,i+2] == letter_vector[2]))
		count[3,i] = length(which(nucleotides[,i+2] == letter_vector[3]))
		count[4,i] = length(which(nucleotides[,i+2] == letter_vector[4]))
	}
	rownames(count) <- c("A", "C", "G", "T")
	return(count)
}
