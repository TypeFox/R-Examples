print.SIMoNeNet <-
function (x, nlimit=20, ...)
{
	cat("Object of class SIMoNeNet (package stringgaussnet)\n\n")
	Nodes <- unique(c(x$Edges$node1,x$Edges$node2))
	cat("Number of nodes:",length(Nodes),"\n")
	Interactions <- unique(paste(x$Edges$node1,x$Edges$node2,sep="."))
	cat("Number of interactions:",length(Interactions),"\n\n")
	cat("Edges preview:\n")
	print(head(x$Edges,nlimit))
	cat("\nDEGenes preview:\n")
	print(head(x$DEGenes,nlimit))
	if (!is.null(x$Annotations))
	{
		cat("\nAnnotations preview:\n")
		print(head(x$Annotations,nlimit))
	}
}
