fisher.TradPerm <-
function (genotypeLine, affectionLine, fromCol, naString, sep, repeatNum = 1000) 
{
    # if (repeatNum<0 || repeatNum != round(repeatNum) ) {
        # stop("'repeatNum' must be a positive integer.")
    # }
    colNum = ncol(genotypeLine)
    genotypeLine = genotypeLine[, fromCol:colNum,drop=FALSE]
    affectionLine = affectionLine[, fromCol:colNum,drop=FALSE]
    # Calculate genotype frequecy for case and control samples
    obsCount = as.matrix(table(affectionLine, genotypeLine))
	colName=colnames(obsCount)
	# delete NA data
    del=0
	if(any(colName==naString)){
	      naIndex=which(colName==naString)
		  obsCount=obsCount[,-naIndex,drop=FALSE]
		  del=1
	}
	
	# test for true data
    fisherP = matrix(0, nrow = 1, ncol = repeatNum)
    obs = fisher.test(obsCount)$p.value
    
	# permute data,then test 
    for (i in 1:repeatNum ) {
	    randIndex=sample(colNum-fromCol+1)
		randData=genotypeLine[,randIndex]
        # Calculate genotype frequecy for case and control samples
        expCount = as.matrix(table(affectionLine, randData))
	    # delete NA data
		if(del==1){
		    expCount=expCount[,-naIndex,drop=FALSE]
		}
        fisherP[1,i] = fisher.test(expCount)$p.value
    }

	# return result
    minExp = fisherP[fisherP < obs]
    pValue = length(minExp)/repeatNum
    list(pValue = pValue, obsP = obs, permP = fisherP)
}
