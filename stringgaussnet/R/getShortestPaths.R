getShortestPaths <-
function (Network,SelectedGenes=0)
{
	optionValue<-getOption("stringsAsFactors")
	options(stringsAsFactors=F)
	if (SelectedGenes==0)
	{
		SelectedGenes<-rownames(Network$DEGenes)
	}
	CombinedScores<-Network$Edges[which(Network$Edges$Interaction=="combined_score"),c("node1","node2","Score")]
	BetweenNodesDistances=convertToDistGraph(CombinedScores) 
	ShortPathsNetwork=findShortestPathways(BetweenNodesDistances,SelectedGenes) 
	FinalNetwork<-ShortPathSTRINGNet(ShortPathsNetwork,Network$DEGenes,Network$Annotations)
	options(stringsAsFactors=optionValue)
	return(FinalNetwork)
}
