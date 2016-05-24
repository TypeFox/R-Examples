mergeFactorEdges <-
function (Networks)
{
	Levels<-names(Networks)
	for (Level in Levels)
	{
		Edges<-Networks[[Level]]$Network$Edges[,c("node1","node2","Rho","P.Value")]
		Edges$Level<-Level
		if(Level==Levels[1]){MergedEdges<-Edges}
		if(Level!=Levels[1]){MergedEdges<-rbind(MergedEdges,Edges)}
	}
	return(MergedEdges)
}
