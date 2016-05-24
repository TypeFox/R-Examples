DHDdist <- function(seq1, l1, seq2, l2, maxcost, sm, norm) {
	seqcost <- 0
	for (l in 1:l1) {
		seqcost <- seqcost + sm[seq1[l]+1,seq2[l]+1,l]
	}
	return(normdist(seqcost, maxcost, l1, l2, norm))
}