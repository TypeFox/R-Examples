getCountsInWindow <-
function(events, startE, endE, windowSize=10000, sorted=FALSE) {
	histBreaks = seq(from=startE, to=endE, by=windowSize)
	if(tail(histBreaks, 1) != endE)
		histBreaks = c(histBreaks, endE)
	countsWindow = hist(events, breaks=histBreaks, right=FALSE, plot=FALSE)$counts
	return(countsWindow)
}

