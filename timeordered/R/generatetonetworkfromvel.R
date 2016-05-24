generatetonetworkfromvel <-
function(vel)
{
	g <- graph.data.frame(vel$edgelist, directed=TRUE, vertices=vel$vertexlist)
	return(g)
}

