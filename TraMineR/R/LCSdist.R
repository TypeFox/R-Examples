## LCS distance between 2 sequences

LCSdist <- function(seq1,l1,seq2,l2,norm) {
	result <- .C(TMR_cLCS, as.integer(seq1), as.integer(seq2), as.double(c(l1,l2)), result = as.integer(0), 
		NAOK=TRUE)$result
	dist <- l1+l2-2*result
	return(normdist(dist, l1+l2, l1, l2, norm))

}
