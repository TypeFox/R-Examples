## ====================================
## Calculates the Length of the Longest
## Common Prefix
## ====================================

seqLLCP <- function(seq1,seq2) {

	if (!inherits(seq1,"stslist") | !inherits(seq2,"stslist")) {
		stop(" [!] sequences must be sequence objects")
	}
	a1 <- alphabet(seq1)
	a2 <- alphabet(seq2)
	if(length(a1)!=length(a2) || any(a1!=a2)){
		stop(" [!] The alphabet of both sequences have to be same.")
	}

	l1 <- seqlength(seq1)
	l2 <- seqlength(seq2)

	result <- .C(TMR_cLCP, as.integer(seq1), as.integer(seq2), as.double(c(l1,l2)), result = as.integer(0))$result

	return(result)
}







