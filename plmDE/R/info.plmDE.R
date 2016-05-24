info.plmDE <- function(Object, ...) {
	cat("Number of samples:", length(Object$sampleInfo[,1]))
	cat("\nNumber of variables regarding these Samples:", length(Object$sampleInfo[1,]) - 2)
	cat("\nNumber of gene expression measurements for each sample:", length(Object$expressionValues[,1]))
	cat("\nSample subgroups:", toString(unique(Object$sampleInfo[,2])))
}