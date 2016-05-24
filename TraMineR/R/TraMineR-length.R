
TraMineR.length <- function(seq, void) {

	totlength <- length(seq)
	nvoid <- sum(seq==void, na.rm=TRUE)
	l <- totlength - nvoid

	return(l)

}
