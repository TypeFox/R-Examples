## ========================================================
## Creates a vector of 0 and 1 used for plotting a sequence
## ========================================================

seqgbar <- function(seq, statl, seql) {

	nbstat <- length(statl)
	
	gbar <- vector("integer", seql*nbstat)

	for (j in 1:seql)
		gbar[((j-1)*nbstat)+which(statl==seq[j])] <- 1

	return(gbar)
}
	
