convertAliases <-
function (FormatInteractions)
{
	SubFormatInteractions=FormatInteractions
	Nodes=unique(c(FormatInteractions$node1,FormatInteractions$node2)) 
	for (Node in Nodes)
	{
		if (requireNamespace("org.Hs.eg.db",quietly=TRUE)) {TempObject<-org.Hs.eg.db::org.Hs.eg} else {stop("org.Hs.eg.db package must be installed to use this function")}
		if (requireNamespace("limma",quietly=TRUE)) {SymboleHGNC=limma::alias2Symbol(Node) } else {stop("limma package must be installed to use this function")}
		if (length(SymboleHGNC)>0) 
		{
			SubFormatInteractions$node1[which(SubFormatInteractions$node1==Node)]=SymboleHGNC[1]
			SubFormatInteractions$node2[which(SubFormatInteractions$node2==Node)]=SymboleHGNC[1]
		}
	}
	return(SubFormatInteractions)
}
