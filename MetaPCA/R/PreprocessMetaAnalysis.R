PreprocessMetaAnalysis <- function(DList, cutRatioByMean=.4, cutRatioByVar=.4, doImpute=FALSE, na.rm.pct=.1, na.rm.pct.each=.5, verbose=FALSE) {	
	studies <- names(DList)
	
	DList <- foreach(dat=iter(DList)) %do% {
		dat[rowSums(is.na(dat))<=ncol(dat)*na.rm.pct.each,] #remove if more than some pct is missing in single study
	}
	
	.genes <- foreach(dat=iter(DList)) %do% {rownames(dat)}
	intersect.rec <- function(.list, ...){
		if(length(.list)==1) return(.list[[1]])
		Recall(c(list(intersect(.list[[1]], .list[[2]], ...)), .list[-(1:2)]), ...)
	}
	.genes <- na.omit(intersect.rec(.genes))
	
	printLog(paste("# of matched genes :", length(.genes)), verbose)
	
	DD <- foreach(dat=iter(DList), .combine=cbind) %do% {
		dat[match(.genes,rownames(dat)),]
	}
	
	DD <- DD[rowSums(is.na(DD))<ncol(DD)*na.rm.pct,]
	
	n <- c(0, cumsum(sapply(DList, ncol)))
	for(i in 1:(length(n)-1)) {
		DList[[i]] <- DD[,(n[i]+1):n[i+1]]
	}
	
	printLog(paste("# of genes after na.rm threshold:", nrow(DList[[1]])), verbose)
	
	if(doImpute) {
		requireAll("impute")
		for(i in 1:length(studies)) {
			if(any(is.na(DList[[i]])))
				DList[[i]] <- impute.knn(DList[[i]])$data
		}
	}
	
	numLeft <- floor(nrow(DList[[1]])*(1-cutRatioByMean)) 
	
	#rank sum of row means
	DList.rankM <- as.matrix(foreach(dat=iter(DList), .combine=cbind) %do% {
				rank(rowMeans(dat))
			})
	DList <- foreach(dat=iter(DList)) %do% {
		dat[order(rowSums(DList.rankM),decreasing=TRUE)[1:numLeft],]
	} 
	
	numLeft <- floor(numLeft*(1-cutRatioByVar))
	
	#rank sum of row var
	DList.rankV <- as.matrix(foreach(dat=iter(DList), .combine=cbind) %do% {
				rank(rowVars(dat))
			})
	DList <- foreach(dat=iter(DList)) %do% {
		dat[order(rowSums(DList.rankV),decreasing=TRUE)[1:numLeft],]
	} 
	
	names(DList) <- studies
	return(DList)
}
