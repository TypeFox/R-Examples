getStatistics <- function(cpm) {
    if (cpm@n < cpm@startup) {return(rep(0,cpm@n))}
    return(getTestStatistics(cpm)$Ds)
}
	

