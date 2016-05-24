rank.partitionedRAD <-
function(radMat, minTrees = 10, min.overall.diff.lnL = 5, 
		 threshold.lnL = 2, discardDoubleCounts = TRUE) {
  radMat.preLnLDiff <- radMat.preMinTrees <- NA
  if(minTrees > 1) {
    nTrees <- apply(radMat, 1, function(x) length(unique(x)))
	radMat.preMinTrees <- radMat
	radMat <- radMat[nTrees >= minTrees, ]
	}
  if(min.overall.diff.lnL > 0) {
    lnL.diff <- apply(radMat, 1, function(x) abs(diff(range(x))))
	radMat.preLnLDiff <- radMat
	radMat <- radMat[(lnL.diff >= min.overall.diff.lnL), ]
	}
  bestMat <- t(apply(radMat, 1, function(x) abs(x - max(x)) <= threshold.lnL))
  worstMat <- t(apply(radMat, 1, function(x) abs(x - min(x)) <= threshold.lnL))
  doubleCountMat <- bestMat & worstMat
  if(discardDoubleCounts) {
    doubleCounts <- apply(doubleCountMat, 1, sum) > 0
	bestMat <- bestMat[!doubleCounts, ]
	worstMat <- worstMat[!doubleCounts, ]
	radMat.preDoubleCounts <- radMat
	radMat <- radMat[!doubleCounts, ]
	}
  out <- list(bestMat = bestMat, worstMat = worstMat, doubleCountMat = doubleCountMat, 
			  params = c(minTrees, min.overall.diff.lnL, threshold.lnL, discardDoubleCounts))
  class(out) <- 'rankedPartitionedRAD'
  out
  }
