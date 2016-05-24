head.DEresults <- function(x, ...) {
	cat("Testing Reults: \n")
	print(head(x$allgenes, ...))
	cat("\n",length(x$DEgenes[,1]), "genes show significant evidence of differential expression. \n")
}
