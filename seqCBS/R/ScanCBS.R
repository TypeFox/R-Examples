ScanCBS <-
function(cases, controls, statistic="binomial", grid.size="auto", takeN=5, maxNCut=100, minStat=0, alpha=0.05, verbose=FALSE, timing=TRUE) {
	## statistic should be either "binomial" (use exact binomial distribution)
	## "normal" (use normal approximation for all windows), or "rabinowitz" (use statistic in Rabinowitz 1994)
	## grid.size should be either "auto" (default) or a vector of integers such as c(10000,1000,100,10)
	timeCBSTotal = proc.time()[3]
	timeCBSPreProcess = proc.time()[3]
	
	nCas = length(cases)
	nCon = length(controls)
	nTotal = nCas + nCon
	p = nCas/nTotal
	lik0 = nCas*log(p) + nCon*log(1-p)
	if(verbose)	{
		print(paste("nCas:", nCas, " nCon:", nCon, " nTotal:", nTotal))
	}
	combCasCon = CombineCaseControlC(cases, controls)
	combL = combCasCon$combL
	combX = combCasCon$combX
	combZ = combCasCon$combZ
	nL = length(combL)
	mBIC = numeric(0)
	curBIC=0
	if(grid.size=="auto") {
		grid.size = getAutoGridSize(nL)
	}
	else if(nL/max(grid.size) > 100) {
		print("Grid is too fine; may result in slow computation")
	}
	grid.size = sort(grid.size, decreasing=TRUE)
	nGridSize = length(grid.size)
	
	min.grid.size.T2 = 2*min(grid.size)
	if(verbose) {
		print(grid.size)
	}
	
	timeCBSPreProcess = proc.time()[3] - timeCBSPreProcess
	
	# tauHat: the vector of change points
	# statHat: matrix of statistics for the change points
	# maxS: matrix with 5 columns, temporary bag of change point candidates
	#       column 1 and 2 are the start and end index of change point
	#       column 3 is the computed statistic with this change point
	#       column 4 and 5 are the left and right bounds of the region for the change points
	
	timeIGSTotal = proc.time()[3]
	timeIGSBreakDown = matrix(0, nrow=nGridSize, ncol=3)
	colnames(timeIGSBreakDown) = c("newComp", "refineComp", "misc")
	rownames(timeIGSBreakDown) = grid.size
	
	tauHat = c(1, nL)
	maxSRes = ScanIterateGrid(combX, combZ, combL, statistic, grid.size, nGridSize, timeIGSBreakDown, takeN, verbose, timing)
	maxS = matrix(c(maxSRes$cptsRet, 1, nL), ncol=5)
	timeIGSBreakDown = maxSRes$timeIGSBreakDown
	statHat = matrix(nrow=0, ncol=6)
	cutStat = minStat + 1
	nCut = 1
	
	timeIGSTotal = proc.time()[3] - timeIGSTotal
	timeCBSBreakDown = c()
	
	while((nCut < maxNCut) && (abs(cutStat) > minStat)) {
		timeCBSnew = proc.time()[3]
		if(verbose) {
			print(dim(maxS))
			print(maxS)
			print(maxS[(max(c(maxS[,2]-maxS[,1], maxS[,1]-maxS[,4], maxS[,5]-maxS[,2])) >= min.grid.size.T2),])
		}
		maxS = matrix(maxS[(max(c(maxS[,2]-maxS[,1], maxS[,1]-maxS[,4], maxS[,5]-maxS[,2])) >= min.grid.size.T2),], ncol=5)
		if(verbose) {
			print(maxS)
		}
		if(length(maxS)==0) {
			break
		}
		maxS.ind = which.max(abs(maxS[,3]))
		maxS.cut = maxS[maxS.ind,]
		mBIC = ScanBIC(combX, combZ, sort(c(tauHat, maxS.cut[1:2])), lik0, nTotal)
		if(mBIC < curBIC && nCut !=1) {
			break
		}
		tauHat = c(tauHat, maxS.cut[1:2])
		statHat = rbind(statHat, c(maxS.cut, mBIC))

		cutStat = maxS.cut[3]
		maxS.cut.L = maxS.cut[4]
		maxS.cut.C1 = maxS.cut[1]
		maxS.cut.C2 = maxS.cut[2]
		maxS.cut.R = maxS.cut[5]
		
		grid.size.L = grid.size[2*grid.size <= maxS.cut.C1 - maxS.cut.L]
		grid.size.M = grid.size[2*grid.size <= maxS.cut.C2 - maxS.cut.C1]
		grid.size.R = grid.size[2*grid.size <= maxS.cut.R - maxS.cut.C2]
		
		timeIGSNew = proc.time()[3]
		
		if(length(grid.size.L) > 0) {
			maxS.LRes = ScanIterateGrid(combX[maxS.cut.L:maxS.cut.C1], combZ[maxS.cut.L:maxS.cut.C1], combL[maxS.cut.L:maxS.cut.C1], statistic, grid.size.L, nGridSize, timeIGSBreakDown, takeN, verbose, timing)
			maxS.L = maxS.LRes$cptsRet
			timeIGSBreakDown = maxS.LRes$timeIGSBreakDown
			maxS.L[c(1,2)] = maxS.L[c(1,2)] + maxS.cut.L -1
			maxS.L = c(maxS.L, maxS.cut.L, maxS.cut.C1)
			maxS = rbind(maxS, maxS.L)
		}
		if(length(grid.size.M) > 0) {
			maxS.MRes = ScanIterateGrid(combX[maxS.cut.C1:maxS.cut.C2], combZ[maxS.cut.C1:maxS.cut.C2], combL[maxS.cut.C1:maxS.cut.C2], statistic, grid.size.M, nGridSize, timeIGSBreakDown, takeN, verbose, timing)
			maxS.M = maxS.MRes$cptsRet
			timeIGSBreakDown = maxS.MRes$timeIGSBreakDown
			maxS.M[c(1,2)] = maxS.M[c(1,2)] + maxS.cut.C1 -1
			maxS.M = c(maxS.M, maxS.cut.C1, maxS.cut.C2)
			maxS = rbind(maxS, maxS.M)
		}
		if(length(grid.size.R) > 0) {
			maxS.RRes = ScanIterateGrid(combX[maxS.cut.C2:maxS.cut.R], combZ[maxS.cut.C2:maxS.cut.R], combL[maxS.cut.C2:maxS.cut.R], statistic, grid.size.R, nGridSize, timeIGSBreakDown, takeN, verbose, timing)
			maxS.R = maxS.RRes$cptsRet
			timeIGSBreakDown = maxS.RRes$timeIGSBreakDown
			maxS.R[c(1,2)] = maxS.R[c(1,2)] + maxS.cut.C2 -1
			maxS.R = c(maxS.R, maxS.cut.C2, maxS.cut.R)
			maxS = rbind(maxS, maxS.R)
		}
		if(verbose) {
			print("maxS after cut")
			print(maxS)
		}
		maxS = matrix(maxS[-maxS.ind,], ncol=5)
		nCut = nCut+1
		curBIC = mBIC
		
		timeIGSTotal = timeIGSTotal + proc.time()[3] - timeIGSNew
		timeCBSBreakDown = c(timeCBSBreakDown, proc.time()[3]-timeCBSnew)
	}
	tauHatInd = sort(tauHat)
	tauHat = combL[tauHatInd]
	
	relCNRes = relCNComp(combX, combZ, tauHatInd, p, alpha)
	relCN = relCNRes$relCN
	relGainLoss = relCNRes$relGainLoss
		
	statHat = cbind(matrix(statHat[,1:2],nrow=nrow(statHat)), statHat)
	statHat[,1] = combL[statHat[,1]]
	statHat[,2] = combL[statHat[,2]]
	statHat[,6] = combL[statHat[,6]]
	statHat[,7] = combL[statHat[,7]]
	colnames(statHat) = c("cptL", "cptR", "cptLReadInd","cptRReadInd", "stat", "parentL", "parentR", "BIC")
	
	timeCBSTotal = proc.time()[3] - timeCBSTotal
	
	if(timing == TRUE) {
		timingRes = list(timeCBSTotal=timeCBSTotal, timeCBSPreProcess=timeCBSPreProcess, timeCBSBreakDown=timeCBSBreakDown, timeIGSTotal=timeIGSTotal, timeIGSBreakDown=timeIGSBreakDown)
		return(list(tauHat=tauHat, statHat = statHat, relCN = relCN, relGainLoss=relGainLoss, timingRes=timingRes))
	}
	else {
		return(list(tauHat=tauHat, statHat = statHat, relCN = relCN, relGainLoss=relGainLoss))
	}
}

