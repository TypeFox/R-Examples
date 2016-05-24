nhppSimulate <-
function(smoothRates) {
	## smoothRate should be a list containing x, y
	## where x: the mid points of the smoothing windows
	## and   y: the smooth number of events in this window
	grid.mid = smoothRates$x
	gridSize = grid.mid[2] - grid.mid[1]
	grid.start = grid.mid - gridSize/2
	smoothRatesY = smoothRates$y
	smoothMaxY = unlist(sapply(1:length(grid.start), function(i) {return(max(smoothRatesY[max(i-1,1)], smoothRatesY[i], smoothRatesY[min(i+1,length(grid.start))]))}))
	hppEvents = unlist(sapply(1:length(grid.start), function(i) {return(sort(hppSimulate(smoothMaxY[i], gridSize)) + grid.start[i])}))
	hppEvents = floor(hppEvents)
	hppEventsRateI = findInterval(hppEvents, grid.start)
	nEvents = length(hppEvents)
	hppEventsIL = findInterval(hppEvents, grid.mid)
	hppEventsIL[hppEventsIL==0] = 1
	hppEventsIR = findInterval(-hppEvents[nEvents:1], -grid.mid[length(grid.mid):1])
	hppEventsIR[hppEventsIR==0] = 1
	hppEventsIR = length(grid.mid)-hppEventsIR+1
	hppEventsIR = hppEventsIR[nEvents:1]
	hppEventsRates = rep(0, nEvents)
	hppEventsProbs = rep(0, nEvents)
	for(i in 1:nEvents) {
		if(hppEventsIL[i] == hppEventsIR[i])
			hppEventsRates[i] = smoothRatesY[hppEventsIL[i]]
		else {
			distL = abs(hppEvents[i]-grid.mid[hppEventsIL[i]])
			distR = abs(hppEvents[i]-grid.mid[hppEventsIR[i]])
			hppEventsRates[i] = (smoothRatesY[hppEventsIL[i]]*distR + smoothRatesY[hppEventsIR[i]]*distL)/(distL+distR)
		}
		temp = hppEventsRates[i]/smoothMaxY[hppEventsRateI[i]]
		if(length(temp)>0) {
			hppEventsProbs[i] = temp
		}
		else {
			hppEventsProbs[i] = 0
		}
	}
	nhppUnifs = runif(nEvents)
	nhppEvents = hppEvents[nhppUnifs < hppEventsProbs]
	return(nhppEvents)
}

