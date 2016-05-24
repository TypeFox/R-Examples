####
# Get statistical information about some R Graph:
####
summary_graph = function(graph_object, graphical_package="igraph") {

	if (graphical_package=="graph") {
		the_graph = graph_object
		# 1. Degree Distribution:
		degree_distribution = graph::degree(graph_object)

		# 2. Nodes per Connected Component:
		print("The analysis of connected components is only possible for undirected graphs; directed graphs will be converted to undirected by RBGL.")
		modules = connectedComp(graph_object)

		nodes_per_module=NULL; j=0
		for (i in modules) {
			j = j + 1
			nodes_per_module[j] = length(i)
		}

		# Subgraphs per connected component:
		graph_per_module = list(); j = 0
		for (i in modules) {
			j = j + 1
			graph_per_module[[j]] <- subGraph(i, graph_object)
		}
	} else {
		the_graph = summary(graph_object)
		degree_distribution = igraph::degree(graph_object)

		modules = clusters(graph_object)
		nodes_per_module=NULL; j=0
		for (i in modules) {
			j = j + 1
			nodes_per_module[j] = length(i)
		}

		graph_per_module = list(); j = 0
		for (i in modules) {
			j = j + 1
			graph_per_module[[j]] <- subgraph(graph_object, i)
		}

		plot(graph_object, layout=layout.fruchterman.reingold, vertex.size=3, vertex.color="green", frame=TRUE, main="Graph (igraph)", vertex.label=NA)
	}

	# 3. Output:
	output = list(nodes_and_edges=the_graph, degree_distribution=degree_distribution, nodes_per_connected_component=nodes_per_module, graph_per_module=graph_per_module)

}
