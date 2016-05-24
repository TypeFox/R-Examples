ScanStatRefineComp <-
function(combZCumSum, combXCumSum, combZPoint, combXPoint, p, nTotal, grid.cur, grid.LR, max.win, statistic) {
	if(statistic=="rabinowitz") {
		SijFactorR = p*(1-p)*(1-1/(nTotal-1))
		#refineRes = ScanStatRefineCompRabin(combZCumSum, combXCumSum, combZPoint, combXPoint, SijFactorR, p, nTotal, grid.cur, grid.L, grid.R, max.win)
		refineRes = .Call("ScanStatRefineCompRabinC", combZCumSum, combXCumSum, combZPoint, combXPoint, SijFactorR, p, nTotal, grid.cur, grid.LR, max.win, PACKAGE="seqCBS")
	}
	else if(statistic=="normal") {
		SijFactorN = p*(1-p)
		refineRes = .Call("ScanStatRefineCompNormalC", combZCumSum, combXCumSum, combZPoint, combXPoint, SijFactorN, p, nTotal, grid.cur, grid.LR, max.win, PACKAGE="seqCBS")
	}
	else if(statistic=="binomial") {
		refineRes = .Call("ScanStatRefineCompBinomC", combZCumSum, combXCumSum, combZPoint, combXPoint, p, nTotal, grid.cur, grid.LR, max.win, PACKAGE="seqCBS")
	}
	else {
		print("Name of Statistic is not Recognized, Rabinowitz used")
		SijFactorR = p*(1-p)*(1-1/(nTotal-1))
		refineRes = .Call("ScanStatRefineCompRabinC", combZCumSum, combXCumSum, combZPoint, combXPoint, SijFactorR, p, nTotal, grid.cur, grid.LR, max.win, PACKAGE="seqCBS")
	}
	return(refineRes)
}

