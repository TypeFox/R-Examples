## Levenshtein distance between seq1 and seq2

levenshtein <- function(seq1, l1, seq2, l2, indel, sm, alphsize, norm) {
	dist <- .C(TMR_cLEVEN, as.integer(seq1), as.integer(seq2), 
		as.double(c(l1,l2,indel,alphsize)), as.double(sm), result = as.double(0),
		NAOK=TRUE)$result
	
	if (norm==4) {
		maxpossiblecost <- (l1+l2)*indel
	}
	else {
		maxpossiblecost <- abs(l1-l2)*indel+min(max(sm), 2*indel)*min(l1,l2)
	}
  return(normdist(dist,maxpossiblecost , l1,l2,norm))
}
