convertToDistGraph <-
function(Network)
{
	AllNodes=unique(c(Network$node1,Network$node2))
	BetweenNodesDistances=matrix(0,ncol=length(AllNodes),nrow=length(AllNodes)) 
	rownames(BetweenNodesDistances)=colnames(BetweenNodesDistances)=AllNodes 
	for (i in 1:nrow(Network))
	{
		BetweenNodesDistances[Network$node1[i],Network$node2[i]]=max(Network$Score)+1-Network$Score[i] 
	}
	return(BetweenNodesDistances)
}
