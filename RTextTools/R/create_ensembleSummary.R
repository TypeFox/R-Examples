create_ensembleSummary <- function(document_summary) {
	label <- function(x) {
		return(paste("n >=",x))
	}
	
	summary <- c()
	for (threshold in 1:max(document_summary$CONSENSUS_AGREE)) {
		algorithms <- document_summary[document_summary$CONSENSUS_AGREE>=threshold,]
		agreement <- round(dim(algorithms)[1]/dim(document_summary)[1],2)
		recall <- round(recall_accuracy(algorithms$MANUAL_CODE,algorithms$CONSENSUS_CODE),2)
		
		summary <- append(summary,c(agreement,recall))
	}

	summary <- matrix(summary,byrow=TRUE,ncol=2)
	colnames(summary) <- c("n-ENSEMBLE COVERAGE","n-ENSEMBLE RECALL")
	rownames(summary) <- c(1:max(document_summary$CONSENSUS_AGREE))
	rownames(summary) <- sapply(rownames(summary),label)

	return(summary)
}