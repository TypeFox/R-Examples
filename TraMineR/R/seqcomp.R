## =======================
## Comparing two sequences
## =======================

seqcomp <- function(x, y) {

	lx <- seqlength(x)		
	ly <- seqlength(y)

	if (lx!=ly) 
		return(FALSE)
	else	if (sum(x[1:lx]==y[1:ly])==lx) 
		return(TRUE)
	else 
		return(FALSE)
}	
