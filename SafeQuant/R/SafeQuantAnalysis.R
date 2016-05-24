# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

### s3 class for storing p-values, ratios and CVs
# created from Expression set (name of controlCondition is specified), adjust. 
# normalization
# replace missing values
# c("global","naRep","rt","pairwise","quantile")
#' safeQunat s3 class
#' @param eset ExpressionSet
#' @param method c("global","naRep","pairwise")
#' @export
safeQuantAnalysis <- function(eset=eset, method=c("global","naRep","pairwise")){
	
	out <- list()
	class(out) <- "safeQuantAnalysis"
	
	### what about fdr, decoy, massTol filters etc..
	#eset <- addIdQvalues(eset)
	
	#normalizationFactors <- NA
	### normalize
	eset <- sqNormalize(eset, method=method)
	
	baselineIntensity <- NA
	# replace missing values
	if("naRep" %in% method ){
		baselineIntensity <- getBaselineIntensity(as.vector(unlist(exprs(eset)[,1])),promille=5)
		### add pseudo (baseline) intensity
		exprs(eset)[is.na(exprs(eset)) | (exprs(eset) < 0)  ] <- 0 
		exprs(eset) <- exprs(eset) + baselineIntensity
		
	}
	
	out$eset <- eset # should the ExpressionSet be stored?
	out$cv <- getAllCV(eset)
	out$ratio <- getRatios(eset,log2=T)
	
	### do not perform stat test for filtered out features
	sel <- !fData(eset)$isFiltered
	
	out$pValue <- out$ratio
	out$pValue[rownames(out$pValue),] <- NA
	
	### we need at least two conditions
	if(length(unique(pData(eset)$condition)) > 1){
		out$pValue[rownames(eset)[sel],] <- getAllEBayes(eset[sel,],adjust=F,method=method)	
	}
	
	out$qValue <- out$ratio
	out$qValue[rownames(out$qValue),] <- NA
	
	### we need at least two runs
	if(length(unique(pData(eset)$condition)) > 1){
		out$qValue[rownames(eset)[sel],] <- getAllEBayes(eset[sel,],adjust=T,method=method)
	}
	out$baselineIntensity <- baselineIntensity
	
	return(out)
	
}


#' @export
.filterSQA <- function(sqa,filter=NA ){
	
	sqa$eset <- sqa$eset[filter,]
	
	tmpNames <- names(sqa$pValue )
	sqa$pValue <- data.frame(sqa$pValue[filter,],row.names=rownames(sqa$eset))
	names(sqa$pValue) <- tmpNames
	
	tmpNames <- names(sqa$qValue )
	sqa$qValue <- data.frame(sqa$qValue[filter,],row.names=rownames(sqa$eset))
	names(sqa$qValue) <- tmpNames
	
	tmpNames <- names(sqa$ratio )
	sqa$ratio <- data.frame(sqa$ratio[filter,],row.names=rownames(sqa$eset))
	names(sqa$ratio) <- tmpNames
	
	tmpNames <- names(sqa$cv )
	sqa$cv <- data.frame(sqa$cv[filter,],row.names=rownames(sqa$eset))
	names(sqa$cv) <- tmpNames
	
	return(sqa)
	
}


#' Export content of safeQuantAnalysis object
#' @param sqa safeQuantAnalysis object
#' @param nbRows Number of rows to export. Features are ordred by increasing minimal p.value
#' @param file file path
#' @export
#' @import limma Biobase
#' @importFrom utils write.table
#' @note  No note
#' @details NA
#' @references NA
#' @seealso  \code{\link{safeQuantAnalysis}}
#' @examples print("No examples")
export <- function(sqa,nbRows=nrow(sqa$pValue), file=NA){
	
	o <- order(apply(sqa$pValue,1,min))
	
	# add pvalue, ratio , cv
	out <- cbind(sqa$pValue,sqa$qValue,sqa$ratio,sqa$cv)
	names(out) <- c(
			paste("pValue",names(sqa$pValue),sep="_")
			,paste("qValue",names(sqa$pValue),sep="_")
			,paste("log2ratio",names(sqa$pValue),sep="_")
			,paste("cv",names(sqa$cv),sep="_")
	
	)
	
	# add feature annotations
	out <- cbind(fData(sqa$eset),exprs(sqa$eset),out)
	
	### export to file
	if(!is.na(file)){
		write.table(out,row.names=F,sep="\t",file=file)
	}else{
		# order features by increasing minimal p-value
		#print(o[1:nbRows])
		print(out[o[1:nbRows],-1])
	}
	
}

#print <- function(x) UseMethod("print")

#' @export
print.safeQuantAnalysis <- function(x, ... ){
	
	#@TODO
	cat("Experimental Design:\n")
	#print(pData(sqa$eset))
	cat("\n")
	
	cat("\nStatistical Analysis:\n")
	print(export(x,nbRows=10))
	if(!is.na(x$baselineIntensity)){
		cat("\nBaseline Intensity:\n")
		cat(x$baselineIntensity,"\n")
	}
	
}

#plot <- function(x) UseMethod("plot")

#' @export
plot.safeQuantAnalysis <- function(x,... ){}