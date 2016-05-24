print.DEGeneExpr <-
function (x, nlimit=20, ...)
{
	cat("Object of class DEGeneExpr (package stringgaussnet)","\n\n")
	cat("Number of samples:",nrow(x$DataExpression),"\n")
	cat("Number of genes:",ncol(x$DataExpression),"\n\n")
	cat("DataExpression preview:\n",sep="")
	print(head(x$DataExpression,nlimit))
	cat("\nDEGenesResults preview:\n",sep="")
	print(head(x$DEGenesResults,nlimit))
}
