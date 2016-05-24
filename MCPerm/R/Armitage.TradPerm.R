Armitage.TradPerm <-
function (genotypeLine, affectionLine, fromCol, naString, sep, repeatNum = 1000) 
{
    # if ( repeatNum<0 || repeatNum != round(repeatNum) ) {
        # stop("'repeatNum' must be a positive integer.")
    # }
	
	# calculate the genotype frequency for case and control samples
    obsStat = genotypeStat(genotypeLine = genotypeLine, affectionLine = affectionLine, fromCol = fromCol, naString = naString, sep = sep)
    obsCount = obsStat$genotypeCount
	
	trendValue = matrix(0, nrow = 1, ncol = repeatNum )
    trendP = matrix(0, nrow = 1, ncol = repeatNum )
	
	# calculate Armitage test for trend in true data
	temp = Armitage(obsCount[[1]],obsCount[[2]],obsCount[[3]],obsCount[[5]],obsCount[[6]],obsCount[[7]])
    obsTrendValue = temp$statistic
    obsTrendP = temp$pValue
	
	N = sum(obsCount)-obsCount[[4]]-obsCount[[8]]
    R = obsCount[[1]]+obsCount[[2]]+obsCount[[3]]
    S = N-R
    n1=obsCount[[2]]+obsCount[[6]]
	n2=obsCount[[3]]+obsCount[[7]]
	temp2 = n1+2*n2
	temp3 = n1+4*n2
	low=S*R*(N*temp3-temp2*temp2)
	
	colNum=ncol(genotypeLine)
	for (i in 1:repeatNum ) {
	    randIndex=sample(colNum-fromCol+1)
	    randData=cbind(genotypeLine[,1:(fromCol-1),drop=FALSE],genotypeLine[,randIndex+fromCol-1,drop=FALSE])
        expStat = genotypeStat(genotypeLine = randData, affectionLine = affectionLine, fromCol = fromCol, naString = naString, sep = sep)
        expCount = expStat$genotypeCount
		temp1=expCount[[2]]+2*expCount[[3]]
		upp=N*(N*temp1-R*temp2)*(N*temp1-R*temp2)
		trendValue[1,i]=upp/low
		trendP[1,i]=pchisq(trendValue[1,i],df=1,lower.tail=FALSE)
	 }
	 
     maxExp = trendValue[trendValue > obsTrendValue]
     pValue = length(maxExp)/repeatNum
     list(pValue = pValue, obsStatistic = obsTrendValue, obsP = obsTrendP, 
	    permStatistic = trendValue, permP = trendP)
}
