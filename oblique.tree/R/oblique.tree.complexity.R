`oblique.tree.complexity` <- function(			#	given
	tree,						#		a) an oblique tree grown using oblique.tree()
	subtree.internal.node.names)			#		b) and a vector of node names of internal nodes
{							#	the complexity of this subtree is calculated
	#calculate the complexity of the tree
	if (dim(tree$frame)[1] > 1) {
		#the tree is not a stump

		#calculate a running sum of how many attributes are used throughout the tree...
		tree.complexity <- 1	#penalty is number of attributes used throughout the tree + 1
		for (internal.nodes.index in which(row.names(tree$frame) %in% subtree.internal.node.names)) {
			#count the number of attributes used at this internal node
#			number.of.node.of.interest <- which(row.names(entire.tree$frame) %in% name.of.node.of.interest)
			if (tree$frame$var[internal.nodes.index] == "") {
				#an oblique split is used at this internal node, count the number of non-zero coefficients
				tree.complexity <- tree.complexity + sum(tree$details[[internal.nodes.index]]$coefficients[-1] != 0)
			} else {
				#an axis-parallel split is used at this internal node
				tree.complexity <- tree.complexity + 1
			}
		}
#print(tree.complexity)	
		return(tree.complexity)
	} else {
		#the tree is a stump
		return(1)				#complexity = size = 1
	}
}

