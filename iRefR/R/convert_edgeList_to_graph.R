####
# Generation of graphical objects from iRefIndex/edgeList tables:
####
convert_edgeList_to_graph = function(edgeList, directionality="undirected", graphical_package="igraph") {
	if (graphical_package == "graph") {
		# 1. Remove polymers and other self-loop edges:
		position_self_loop_edges = which(edgeList[,1]==edgeList[,2])
		nodes_in_loops = unique(edgeList[position_self_loop_edges,1])
		if (length(nodes_in_loops)>0) {
			cat("Polymers and loop edges will be removed. If you want to work with a multigraph (loops included), use igraph as graphical_package.\n")
			edgeList = edgeList[-position_self_loop_edges,]	# get edgeList without loops
		}
		cat("RIGID and edgeType edge attributes are not supported for the chosen \"graph\" package.\n")
		edgeList = unique(edgeList[,c(1,2,5)])

		# 2. Generate graph:
		edges = edgeList[,1:2]
		weights = as.integer(edgeList[,3])
		the_graph = ftM2graphNEL(edges, W=weights, edgemode=directionality)
	} else {
		edgel = as.data.frame(edgeList)
		colnames(edgel) = c("V1", "V2", "rigid", "edgeType", "weight")
		if (directionality=="directed") {
			direction_flag = "TRUE"
		} else {
			direction_flag = "FALSE"
		}
		the_graph = graph.data.frame(edgel, directed=direction_flag)
	}

	output = the_graph
}
