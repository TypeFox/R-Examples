processObservation <- function(cpm,x) {
    cpm@n <- cpm@n + 1
	n <- cpm@n #fix

    #if the threshold sequence is not long enough, then expand it
    if (cpm@n > length(cpm@hs)) {
        cpm@hs <- c(cpm@hs, rep(cpm@hs[length(cpm@hs)],10000))
    }
    
	cpm@windowStatistic <- updateWindowStatistic(cpm,x)

    # must be handled separately since both component CPMS 
    # need to be updated
	if (class(cpm)=='ChangePointModelLepage') {
		cpm <- cpmLepageProcessObservation(cpm,x)
	}
    

	#now try to detect a change
	if (length(cpm@hs)>0 && cpm@n >= cpm@startup) { 
		val <- getTestStatistics(cpm)$val
				
		if (!is.na(val) && val > cpm@hs[cpm@n]) {
			cpm@changeDetected <- TRUE
		} 
	}

	return(cpm)
}