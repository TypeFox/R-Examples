####
# Generation of iRefIndex/edgeList tables from graphical objects:
####
convert_graph_to_edgeList = function(graph, directionality="undirected", graphical_package="igraph") {
	if (graphical_package == "graph") {
		# 1. Extract edges and weights:
		edgeList = edges(graph)
		weightList = edgeWeights(graph)

		# 2. Reformat edges:
		n = 0; m = 0; this_edge=list(); this_weight=list()
		for (i in 1:length(edgeList)) {
			tmp = edgeList[[i]]
			weights = weightList[[i]]
			for (j in tmp) {
				n = n+1
				this_edge[[n]] = c(names(edgeList)[[i]], j)
			}
			for (j in weights) {
				m = m+1
				this_weight[[m]] = j
			}
		}
		edge_table = do.call(rbind, this_edge)
		weight_table = do.call(rbind, this_weight)
		the_table = cbind(edge_table, weight_table)

		# 3. Remove repeated edges:
		total_unique = function(edgeList) {# Remove repeated lines counting AB as equal to BA
			edgeList = unique(edgeList)
			ordered_edgeList = list()
			for (i in 1:dim(edgeList)[1]) {
				ordered_edgeList[[i]] = c(sort(edgeList[i,1:2]), edgeList[i,3])
			}
			edges = do.call(rbind, ordered_edgeList)
			edges = unique(edges)
		}
		if (directionality=="directed") {
			final_table = the_table
		} else {
			final_table = total_unique(the_table)
		}

	} else {
		edges = get.edgelist(graph)
		weights = E(graph)$weight
		final_table = cbind(edges, weights)
	}

	output = final_table
}
