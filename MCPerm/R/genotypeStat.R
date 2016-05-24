genotypeStat <-
function (genotypeLine, affectionLine, fromCol, naString, sep) 
{
    # if (nrow(genotypeLine) != 1) {
        # stop("'genotypeLine' must be a data.frame or matrix with nrow=1.")
    # }
	
    # if (nrow(affectionLine) != 1) {
        # stop("'affectionLine' must be a data.frame or matrix with nrow=1.")
    # }
	
    # if (ncol(genotypeLine) != ncol(affectionLine)) {
        # stop("'genotypeLine' and 'affectionLine' must have the same length.")
    # }
	
    # if (fromCol<0 || fromCol != round(fromCol)) {
        # stop("'fromCol' must be a positive integer.")
    # }

    # if (!is.character(naString)) {
        # stop("'naString' must be a character string.")
    # }
    # if (!is.character(sep)) {
        # stop("'sep' must be a character string.")
    # }
	
    colNum = ncol(genotypeLine)
    genotypeCount = genotypeLine[, fromCol:colNum,drop=FALSE]
    affection = affectionLine[, fromCol:colNum,drop=FALSE]
	# Calculate genotype frequecy for case and control samples
    tempResult = as.matrix(table(affection, genotypeCount))
	colName=colnames(tempResult)
	colNum=ncol(tempResult)
	# delete NA data
	if(any(colName==naString)){
	     if(colNum==1){
		     return(tempResult)
		 }
	     naIndex=which(colName==naString)
		 naResult=tempResult[,naIndex,drop=FALSE]
		 tempResult=tempResult[,-naIndex,drop=FALSE]
		 colName=colName[-naIndex]
		 colNum=ncol(tempResult)
	}else{
	     naResult=matrix(0,nrow=2,ncol=1)
	}
	
	# # deal with "GA==AG" and statistic allele number
	# alleleName=c()
    # delIndex = c()
    # for (i in 1:colNum) {
        # colNameSplit = unlist(strsplit(colName[i], sep))
		# alleleName=c(alleleName,colNameSplit)
        # if (colNameSplit[1] != colNameSplit[2]) {
            # colNameSplit = sort(colNameSplit)
            # colName[i] = paste(colNameSplit[1], colNameSplit[2], sep = sep)
            # colNameIndex = which(colName == colName[i])
            # if (length(colNameIndex) == 2) {
                # tempResult[, colNameIndex[1]] = tempResult[, colNameIndex[1]] + tempResult[, colNameIndex[2]]
                # delIndex = c(delIndex, colNameIndex[2])
            # }
        # }
    # }
	# allele=unique(alleleName)
	# if(length(allele)>2){
	     # stop("allele number is more than 2.")
	# }
    # if (!is.null(delIndex)) {
        # tempResult = tempResult[, -delIndex,drop=FALSE]
		# colNum = ncol(tempResult)
    # }
	
	# deal with genotype<3 
	alleleName=c()
	for(i in 1:colNum){
	    alleleName=c(alleleName,unlist(strsplit(colName[i], sep)))
	}
	
	rowNum = nrow(tempResult)
	if(colNum==2){
	     tempCount=matrix(0,nrow=rowNum,ncol=1)
		 if (alleleName[1] != alleleName[2]){
			 tempResult=cbind(tempCount,tempResult)
		 }else{
			 if (alleleName[3] != alleleName[4]){
				 tempResult=cbind(tempResult,tempCount)
			 }else{
                 tempResult=cbind(tempResult[,1,drop=FALSE],tempCount,tempResult[,2,drop=FALSE])
			 }
		 }
	}
	if(colNum==1){
	     tempCount=matrix(0,nrow=rowNum,ncol=1)
		 if (alleleName[1] != alleleName[2]){
			 tempResult=cbind(tempCount,tempResult,tempCount)
		 }else{
		     tempResult=cbind(tempResult,tempCount,tempCount)
		 }
	}
	
	# calculate allele frequecy for case and control
	alleleResult=matrix(0,nrow=rowNum,ncol=2)
	for(r in 1:rowNum){
	     alleleResult[r,1]=2*tempResult[r,1]+tempResult[r,2]
		 alleleResult[r,2]=2*tempResult[r,3]+tempResult[r,2]
	}
	
	# return result
	alleleCount=c('case_allele1'=alleleResult[2,1], 'case_allele2'=alleleResult[2,2],
	              'control_allele1'=alleleResult[1,1],'control_allele2'=alleleResult[1,2])
    genotypeCount=c('case_11'=tempResult[2,1],'case_12'=tempResult[2,2],'case_22'=tempResult[2,3],'case_NA'=naResult[2,1],
	     'control_11'=tempResult[1,1],'control_12'=tempResult[1,2],'control_22'=tempResult[1,3],'control_NA'=naResult[1,1])
	return(list(alleleCount = alleleCount, genotypeCount = genotypeCount))
}
