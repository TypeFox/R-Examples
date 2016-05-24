## ============================
## Number of matching positions
## ============================

seqmpos <- function(seq1, seq2, with.missing=FALSE) {

	if (!inherits(seq1,"stslist") | !inherits(seq2,"stslist")) 
		stop("sequences must be sequence objects")

	## Defining the positions to compare with logical values
	## void positions are set to FALSE and hence will not be
	## counted in the sum of matching positions
	comp1 <- seq1!=attr(seq1,"void")
	comp2 <- seq2!=attr(seq2,"void")

	if (!with.missing) {
		comp1 <- comp1 & seq1!=attr(seq1,"nr")
		comp2 <- comp2 & seq2!=attr(seq2,"nr")
	}

	mpos <- sum(seq1==seq2 & comp1 & comp2)

	return(mpos)
}
