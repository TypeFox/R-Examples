selectGroup <-
function(design, ydata, varNames, nvarGroup, typeRF=ifelse(is.factor(ydata), "classif", "reg"), verbose=TRUE, ntree=500, ...){

	if(length(unique(ydata))<=2 & !is.factor(ydata))
		ydata <- as.factor(ydata)

	typeError <- ifelse(typeRF=="classif", "err.rate", "mse")

	if(missing(varNames) & missing(nvarGroup))
		stop("You must supply 'varNames' or 'nvarGroup")

	if(missing(varNames))
		varNames <- names(nvarGroup)

	n <- length(ydata)
	p <- length(varNames)
 	
 	if(missing(nvarGroup))
 		stop("nvarGroup is missing")

 	if(verbose){
 		cat("normalize =", (length(unique(nvarGroup))!=1), "\n")
 		print(nvarGroup)
 	}

 	indexesVarInDesign <- .extractIndexes(nvarGroup)

	if(verbose) 
		cat(ifelse(typeRF=="reg", "Regression", "Classification"), "backward selection.\nSplitting data into a training and a testing set...")

	u <- sample(c("idxTr", "idxTe"), n, replace=TRUE, prob=c(2/3,1/3))
	spl <- split(1:n, u)
	idxTr <- spl$idxTr
	idxTe <- spl$idxTe

	if(verbose) cat("\n")


	error <- ranking <- rankingIndexes <- numeric(p)
	varIndex <- 1:p
	varNames.iter <- varNames
	i <- p
	while(i>=1){

		if((i!=length(varIndex)) | (i!=length(varNames.iter))){
			cat("i =", i, "length(varIndex) =", length(varIndex), "length(varNames.iter) =", length(varNames.iter), "\n")
			stop("length(varIndex)")
		}
		indexes <- unlist(.extractIndexesSub(varNames.iter, indexesVarInDesign))
		if(verbose) cat("Survival indexes :", indexes, "\n")
	
		if(is.list(indexes)){
			survivalIndexes <- unlist(indexes)
		}else{
			survivalIndexes <- as.numeric(indexes)
		}

		forest <- randomForest(	x=as.matrix(design[idxTr,survivalIndexes]), y=ydata[idxTr], xtest=as.matrix(design[idxTe,survivalIndexes]), 
								ytest=ydata[idxTe], keep.forest=TRUE, keep.inbag=TRUE, ntree=ntree) #, ...)
		error[i] <- forest$test[[typeError]][ntree]

		if(i!=1){

			idxGroup <- 1:length(survivalIndexes) - 1

			IMP <- varImpGroup(object=forest, xdata=design[idxTr,survivalIndexes], ngroups=i, nvarGroup=nvarGroup[varNames.iter], idxGroup=idxGroup, groupsNames=varNames.iter, ...)
			if(verbose) print(sort(round(IMP, 3), decreasing=TRUE))
			minIndex <- which.min(IMP)
			ranking[i] <- varNames.iter[minIndex]
			rankingIndexes[i] <- which(varNames==ranking[i])
			if(verbose) cat(varNames.iter[minIndex], "eliminated.", i-1, "remaining groups of variables.\t Error =", round(error[i], 2), "\n\n\n")
		}else{
			ranking[i] <- varNames.iter
			rankingIndexes[i] <- which(varNames==varNames.iter)
			if(verbose) cat(varNames.iter, "eliminated. No remaining groups of variables.\t Error =", round(error[i], 2), "\n\n\n")
		}
		i <- i-1
		varIndex <- varIndex[-minIndex]
		varNames.iter <- varNames.iter[-minIndex]
	}
	if(verbose) cat("Ending...")
	nselected <- which.min(error)
	selection <- ranking[1:nselected]
	selectionIndexes <- rankingIndexes[1:nselected]

	if(verbose) cat("\n", nselected, "selected variables:\n", selection, "\n\n")

	liste <- list(	nselected=nselected, selection=selection, selectionIndexes=selectionIndexes, error=error, ranking=ranking, rankingIndexes=rankingIndexes, 
					typeRF=typeRF)

	class(liste) <- "fRFE"
	return(liste)
}
