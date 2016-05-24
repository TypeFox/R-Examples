nhppRateEstimate <-
function(controls, length.out=floor(length(controls)/20), lowessF=0.1) {
	maxVal = max(controls)
	grid.fix = seq(1, maxVal, length.out=length.out)
	gridSize = grid.fix[2]-grid.fix[1]
	grid.mid = grid.fix + gridSize/2
	controlCountInGrid = getCountsInWindow(controls, 0, maxVal, gridSize, sorted=FALSE)
	controlSmoothRates = lowess(x=grid.mid, y=controlCountInGrid, lowessF)
	controlSmoothRates$y = controlSmoothRates$y/(gridSize)
	return(controlSmoothRates)
}

