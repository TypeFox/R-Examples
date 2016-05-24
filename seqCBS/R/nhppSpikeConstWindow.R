nhppSpikeConstWindow <-
function(smoothRates, nSpike=25, cptLen=5, minGain=1.5, maxGain=10, minLoss=0.01, maxLoss=0.5, pGain=0.6) {
	grid.mid = smoothRates$x
	spikeRate = smoothRates$y
	gridSize = grid.mid[2] - grid.mid[1]
	grid.fix = grid.mid - gridSize/2
	nGrid = length(grid.fix)
	
	tempGrid = seq(1, (nGrid-cptLen), by=(cptLen+1))
	nSpike = min(length(tempGrid), nSpike)
	gridL = sample(tempGrid, nSpike)
	gridR = gridL+cptLen-1
	gridR[gridR>nGrid] = nGrid
	relCN = sample(0:1, nSpike, replace=TRUE, prob=c(1-pGain, pGain))
	relCN[1] = 0
	relCN[2] = 1
	relCN[relCN==0] = runif(sum(relCN==0), min=minLoss, max=maxLoss)
	relCN[relCN==1] = runif(sum(relCN==1), min=minGain, max=maxGain)
	spikeMat = cbind(gridL, gridR, grid.fix[gridL], grid.fix[gridR]+gridSize, relCN)
	colnames(spikeMat) = c("gridL", "gridR", "readL", "readR", "relCN")
	
	for(i in 1:nSpike) {
		spikeRate[gridL[i]:gridR[i]] = spikeRate[gridL[i]:gridR[i]]*relCN[i]
	}
	meanSpikeRate = mean(spikeRate)
	caseRates = list(x=grid.mid, y=spikeRate)
	return(list(spikeMat = spikeMat, caseRates=caseRates))
}

