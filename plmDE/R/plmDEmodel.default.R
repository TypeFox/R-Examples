plmDEmodel.default <- function(genes, expressionValues, sampleInfo, ...) {
	## expressionValues must be a data frame with one column for each sample
	## and one row for each gene.  sampleInfo must be a data frame
	## where the ith entry in the first column matches the ith sample name
	## in expressionValues.  
	expressionValues = as.data.frame(expressionValues)
	sampleInfo = as.data.frame(sampleInfo)
	if (length(genes) != length(expressionValues[,1])) {
		stop("Length of genes list must match the number of expression values in each column.")
	}
	for (i in 1:length(sampleInfo[,1])) {
		if (sampleInfo[i,1] != names(expressionValues)[i]) {
			stop("entries in sampleInfo[i, 1] must correspond with column names of expressionValues")
		}
	}
	dataObject = list(genes = genes, expressionValues = expressionValues, sampleInfo = sampleInfo)
	class(dataObject) = "plmDE"
	return(dataObject)
}