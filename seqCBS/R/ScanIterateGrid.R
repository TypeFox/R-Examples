ScanIterateGrid <-
function(combX, combZ, combL, statistic, grid.size, nGridSize, timeIGSBreakDown, takeN, verbose, timing) {
	## Argument: combX, combZ, name of statistic, grid size, takeN, verbose
	## Return:
	## Do a variable window scan statistic of the case/control poisson processes
	timeIGSnew = proc.time()[3]
	nCas = sum(combZ)
	nTotal = as.numeric(sum(combX))
	nCon = nTotal-nCas
	nL = length(combL)
	p = nCas/nTotal
	if(verbose)	{
		print(paste("nCas:", nCas, " nCon:", nCon, " nTotal:", nTotal, "nL:", nL))
	}
	timeIGSBreakDown[,3] = proc.time()[3] - timeIGSnew + timeIGSBreakDown[,3]

	
	## Initialize the bag of change point calls
	cpts = matrix(nrow=0, ncol=3)
	if(grid.size[length(grid.size)] > 4) {
		grid.size = c(grid.size, 2)
	}
	
	## Variable Windows by Grid Refinement
	for(g in 1:length(grid.size)) {
		timeRow = nGridSize - length(grid.size) + g
		timeIGSnew = proc.time()[3]
		## 1. Construct the current grid and define max.win
		grid.cur = as.numeric(seq(1, nL, grid.size[g]))
		if(grid.cur[length(grid.cur)] != nL)	grid.cur=c(grid.cur, nL+0.0)
		if(g==1){
			max.win=length(grid)
		}
		else{
			max.win = 2*floor(grid.size[g-1]/grid.size[g])
		}
		max.win=as.numeric(max.win)
		combZCumSum = as.numeric(.Call("ScanIGSGridCumSumC", as.numeric(combZ), as.numeric(grid.cur), PACKAGE="seqCBS"))
		combXCumSum = as.numeric(.Call("ScanIGSGridCumSumC", as.numeric(combX), as.numeric(grid.cur), PACKAGE="seqCBS"))
		combZPoint = as.numeric(combZ[grid.cur])
		combXPoint = as.numeric(combX[grid.cur])
		timeIGSBreakDown[timeRow,3] = timeIGSBreakDown[timeRow,3] + proc.time()[3] - timeIGSnew
		if(verbose)	{
			print(paste("PassedCumSum", g))
		}
		
		## 2. Refine existing change-points in bag
		timeIGSnew = proc.time()[3]
		if(nrow(cpts)>0) {
			grid.LR = cpts %/% grid.size[g]
			for(r in 1:nrow(cpts)) {
				#print(c(cpts[r,], grid.LR[r,], grid.size[g]))
				refineRes = ScanStatRefineComp(combZCumSum, combXCumSum, combZPoint, combXPoint, p, nTotal, grid.cur, as.numeric(grid.LR[r,]), max.win, statistic)
				#print(refineRes)
				cpts[r,] = refineRes[which.max(abs(refineRes[,3])),]
			}
		}
		if(verbose)	{
			print(paste("PassedRefineScan", g))
		}
		timeIGSBreakDown[timeRow,2] = timeIGSBreakDown[timeRow,2] + proc.time()[3] - timeIGSnew
		
		## 3. New Scan with the current grid size
		if(grid.size[g] > 4) {
			timeIGSnew = proc.time()[3]
			newRes = ScanStatNewComp(combZCumSum, combXCumSum, combZPoint, combXPoint, p, nTotal, grid.cur, max.win, statistic)
			if(takeN/nrow(newRes) > 0.2) {
				cpts = rbind(cpts, newRes[order(abs(newRes[,3]), decreasing=TRUE)[1:min(c(takeN, nrow(newRes)))],])
			}
			else {
				abs3NewRes = abs(newRes[,3])
				for(k in 1:takeN) {
					maxIndex = which.max(abs3NewRes)
					abs3NewRes = abs3NewRes[-maxIndex]
					cpts = rbind(cpts, newRes[maxIndex,])
					newRes = newRes[-maxIndex,]
				}
			}
			if(verbose)	{
				print(paste("PassedNewScan", g))
			}
			timeIGSBreakDown[timeRow,1] = timeIGSBreakDown[timeRow,1] + proc.time()[3] - timeIGSnew
		}
	}
	cptsRet = cpts[which.max(abs(cpts[,3])),]
	return(list(cptsRet=cptsRet, timeIGSBreakDown=timeIGSBreakDown))
}

