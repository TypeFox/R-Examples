keepOC <-
function(betaTruth,
									 betaHat,
									 seHat,
									 threshold=c(-Inf,Inf))
{
	## 'value' is a B x nDesigns matrix that indicates which replicates to keep for each design
	
	##
	B        <- dim(betaHat)[1]
	nDesigns <- dim(betaHat)[2]
	value    <- matrix(TRUE, nrow=B, ncol=nDesigns)
	
	##
	for(i in 1:nDesigns)
	{
		##
		prob1 <- apply(is.na(betaHat[,i,-1]), 1, sum) > 0
		prob2 <- apply(is.na(seHat[,i,-1]), 1, sum) > 0
		prob3 <- (betaHat[,i,-1] < array(threshold[1], dim=dim(betaHat[,i,-1]))) + (betaHat[,i,-1] > array(threshold[2], dim=dim(betaHat[,i,-1])))
		prob3 <- apply(prob3 > 0, 1, sum) > 0
		prob3[is.na(prob3)] <- TRUE

		##
		value[,i] <- (prob1 == FALSE & prob2 == FALSE & prob3 == FALSE)
	}
	
	##
	return(value)
}
