generatetimeaggregatednetwork <-
function(g, starttime, stoptime)
{
	sg <- induced.subgraph(g, V(g)[V(g)$Time >= starttime & V(g)$Time < stoptime])

	newedgelist <- data.frame(VertexFrom=E(sg)$VertexFrom,VertexTo=E(sg)$VertexTo,stringsAsFactors=FALSE)
	newedgelist <- newedgelist[newedgelist[,1] != newedgelist[,2],]
	newedgelist$EdgePair <- as.numeric(factor(paste(newedgelist$VertexFrom, 	newedgelist$VertexTo)))
	tabcounts <- tabulate(newedgelist$EdgePair)
	newedgelist$Count <- tabcounts[newedgelist$EdgePair]
	newedgelist$weight <- newedgelist$Count
	
	uniqueedgelist <- unique(newedgelist)
	uniqueedgelist <- uniqueedgelist[,c("VertexFrom","VertexTo","weight","Count")]
	allvertices <- as.data.frame(unique(V(g)$Name))
	
	timeaggregatednetwork <- graph.data.frame(d=uniqueedgelist, vertices=allvertices)

	return(timeaggregatednetwork)
}

