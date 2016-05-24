findShortestPathways <-
function(BetweenNodesDistances,SelectedNodes)
{
	ShortPathDistances=matrix(0,ncol=length(SelectedNodes),nrow=length(SelectedNodes)) 
	colnames(ShortPathDistances)=rownames(ShortPathDistances)=SelectedNodes
	ShortPaths=Intermediates=ShortPathDistances 
	ShortPaths[which(ShortPaths==0)]="" 
	for (i in 1:length(SelectedNodes)) 
	{
		OtherIndices=(i+1):length(SelectedNodes) 
		OtherIndices=OtherIndices[which(OtherIndices!=i & OtherIndices<=length(SelectedNodes))] 
		for (j in OtherIndices)
		{
			source=SelectedNodes[i] 
			target=SelectedNodes[j] 
			OtherSelectedNodes=SelectedNodes[which(!(SelectedNodes %in% c(source,target)))]
			SubDistances=BetweenNodesDistances[which(!(rownames(BetweenNodesDistances) %in% OtherSelectedNodes)),which(!(colnames(BetweenNodesDistances) %in% OtherSelectedNodes))] 
			if (requireNamespace("igraph",quietly=TRUE)) {SubGraph=igraph::graph.adjacency(SubDistances,mode="undirected",weighted=T)} else {stop("igraph package must be installed to use this function")}
			if (requireNamespace("igraph",quietly=TRUE)) {Nodes=igraph::get.vertex.attribute(SubGraph,"name")} else {stop("igraph package must be installed to use this function")}
			sourceid=which(Nodes==source)
			targetid=which(Nodes==target)
			if (requireNamespace("igraph",quietly=TRUE)) {Distance=igraph::shortest.paths(SubGraph,sourceid,targetid)} else {stop("igraph package must be installed to use this function")}
			if (is.finite(Distance)) 
			{
				if (requireNamespace("igraph",quietly=TRUE)) {ShortPath=igraph::get.shortest.paths(SubGraph,sourceid,targetid)$vpath[[1]]} else {stop("igraph package must be installed to use this function")}
				ShortPath=ShortPath[which(!(ShortPath %in% c(sourceid,targetid)))]
				ShortPathDistances[source,target]=ShortPathDistances[target,source]=Distance 
				ShortPaths[source,target]=ShortPaths[target,source]=paste(Nodes[ShortPath],collapse=",") 
				Intermediates[source,target]=Intermediates[target,source]=length(Nodes[ShortPath]) 
			}
		}
	}
	return(list(Distances=ShortPathDistances,ShortPaths=ShortPaths,Intermediates=Intermediates))
}
