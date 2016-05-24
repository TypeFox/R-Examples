print.STRINGNet <-
function(x, nlimit=20, ...)
{
	cat("Object of class STRINGNet (package stringgaussnet)\n\n")
	AllNodes <- unique(c(x$Edges$node1,x$Edges$node2))
	InitNodes <- rownames(x$DEGenes)
	AddedNodes <- AllNodes[which(!(AllNodes %in% InitNodes))]
	cat("Total number of nodes:",length(AllNodes),"\n")
	cat("Number of initial nodes:",length(InitNodes),"\n")
	cat("Number of added nodes:",length(AddedNodes),"\n\n")
	AllInteractions <- unique(paste(x$Edges$node1,x$Edges$node2,sep="."))
	InitInteractions <- unique(paste(x$Edges$node1,x$Edges$node2,sep=".")[which(x$Edges$node1 %in% InitNodes & x$Edges$node2 %in% InitNodes)])
	AddedInteractions <- AllInteractions[which(!(AllInteractions %in% InitInteractions))]
	cat("Total number of interactions:",length(AllInteractions),"\n")
	cat("Number of interactions between initial nodes:",length(InitInteractions),"\n")
	cat("Number of interactions with added nodes:",length(AddedInteractions),"\n\n")
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
