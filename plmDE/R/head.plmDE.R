head.plmDE <- function(x, ...) {
	cat(" Genes: \n")
	print(head(x$genes, ...))
	cat("\nInformation about Samples: \n")
	print(head(x$sampleInfo, ...))
	cat("\nExpression Values: \n")
	print(head(x$expressionValues, ...))
}
