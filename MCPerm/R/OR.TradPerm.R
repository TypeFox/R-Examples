OR.TradPerm <-
function (genotypeLine, affectionLine, fromCol, naString, sep, repeatNum = 1000) 
{
    # if ( repeatNum<0 || repeatNum != round(repeatNum) ) {
        # stop("'repeatNum' must be a positive integer.")
    # }
	
	# calculate the genotype frequency for case and control samples
    obsStat = genotypeStat(genotypeLine = genotypeLine, affectionLine = affectionLine, fromCol = fromCol, naString = naString, sep = sep)
    obsCount = obsStat$alleleCount
	
	ORValue = matrix(0, nrow = 1, ncol = repeatNum)
	
	# calculate OR in true data
	temp = OR(obsCount[[1]],obsCount[[2]],obsCount[[3]],obsCount[[4]])
	risk_allele=temp$risk_allele
    obs = temp$OR
	
	colNum=ncol(genotypeLine)
	for (i in 1:repeatNum) {
	    randIndex=sample(colNum-fromCol+1)
	    randData=cbind(genotypeLine[,1:(fromCol-1),drop=FALSE],genotypeLine[,randIndex+fromCol-1,drop=FALSE])
        expStat = genotypeStat(genotypeLine = randData, affectionLine = affectionLine, fromCol = fromCol, naString = naString, sep = sep)
        expCount = expStat$alleleCount
		if(risk_allele==1){
		     ORValue[1,i]=(expCount[[1]]*expCount[[4]])/(expCount[[2]]*expCount[[3]])
		}else{
		     ORValue[1,i]=(expCount[[2]]*expCount[[3]])/(expCount[[1]]*expCount[[4]])
		}
	}
    maxExp = ORValue[ORValue > obs]
    pValue = length(maxExp)/repeatNum
    list(risk_allele=risk_allele,pValue = pValue, obsOR = obs, permOR = ORValue)
}
