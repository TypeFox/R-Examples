FilterEdges.ShortPathSTRINGNet <-
function (x,Threshold,AttributeFilter="Distance", ...)
{
	KeepPos<-which(x$Edges[,AttributeFilter]<=Threshold)
	FilteredNet<-x
	FilteredNet$Edges<-FilteredNet$Edges[KeepPos,]
	Nodes<-unique(c(FilteredNet$Edges$node1,FilteredNet$Edges$node2))
	FilteredNet$DEGenes<-FilteredNet$DEGenes[rownames(FilteredNet$DEGenes) %in% Nodes,]
	if(!is.null(FilteredNet$Annotations)){FilteredNet$Annotations<-FilteredNet$Annotations[rownames(FilteredNet$Annotations) %in% Nodes,]}
	return(FilteredNet)
}
