"rg.remove.na" <- 
function(xx)
{
	# Function to remove NAs from a vector and inform the user of how many.
	#
	nna <- 0
	wna <- is.na(xx)
	if(any(wna)) {
		x <- xx[!wna]
		nna <- sum(wna)
		cat(" ", nna, "NA (missing) values removed from input vector\n")
	}
	else x <- xx
	invisible(list(x = x, nna = nna))
}

