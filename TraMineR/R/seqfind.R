## =================================
## Find the ocurrences of a sequence
## =================================

seqfind <- function(x,y) {

	if (!inherits(x,"stslist")) 
		stop("x is not a sequence object, use 'seqdef' function to create one")

	if (!inherits(y,"stslist")) 
		stop("y is not a sequence object, use 'seqdef' function to create one")
	
	x.void <- attr(x,"void")
	y.void <- attr(y,"void")

	x.conc <- as.vector(seqconc(x, void=x.void))
	y.conc <- as.vector(seqconc(y, void=y.void))

	occ <- NULL

	for (i in 1:length(x.conc))
		occ <- c(occ,which(y.conc==x.conc[i]))

	return(occ)
}
	 
	
	
