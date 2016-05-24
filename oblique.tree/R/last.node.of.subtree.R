`last.node.of.subtree` <- function(			#	given
	tree,						#		a) a tree grown using oblique.tree() or tree()
	subtree.root.name)				#		b) the name of an internal node
{							#	the row number of the last node in this subtree is found
#temp <- subtree.root.name
	repeat {
		right.child.name <- 2 * subtree.root.name + 1		#got to right child
		does.right.child.exist <- row.names(tree$frame) %in% right.child.name
		if (sum(does.right.child.exist) == 1) {
			#right child exists, move current node of interest to that node
			subtree.root.name <- right.child.name
			current.node.row.number <- which(does.right.child.exist)
		} else {
			break
		}
	}
#print(c(temp,row.names(tree$frame)[current.node.row.number]))
	return(current.node.row.number)
}

