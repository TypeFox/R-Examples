selectLevel <-
function(design, ydata, typeRF=ifelse(is.factor(ydata), "classif", "reg"), verbose=TRUE, ntree=500, ...){
	varNamesLevels <- unique(sapply(colnames(design), FUN=function(str) strsplit(str, "_")[[1]][1] ))
	
	wrap <- function(var){
		tmp <- ifelse(substr(var,1,1) == 's', var, paste(var, '_', sep=''))
		return(length(grep(tmp, colnames(design))))
	}
	nvarGroup <- sapply(varNamesLevels, FUN=wrap)
	if(verbose) cat("Group names:", varNamesLevels, "\tNr of variable in each level:", nvarGroup, "\n")

	selectGroup(design=design, ydata=ydata, varNames=varNamesLevels, nvarGroup=nvarGroup, typeRF=typeRF, verbose=verbose, ntree=ntree, ...)
}
