`oblique.tree.prune.nodes` <- function(			#	given
	tree,						#		a) a tree
	list.of.node.names.to.prune)			#		b) and a list containing the names of nodes all nodes in subtrees
{							#	the tree is updated to reflect the pruning of these subtrees to a leaf
#browser()
	new.leaves <- NULL
	for (list.index in 1:length(list.of.node.names.to.prune)) {
		#edit the row of tree$frame that becomes a leaf
		rows.to.update <- row.names(tree$frame) %in% list.of.node.names.to.prune[[list.index]][1]
		new.leaves <- 	c(	new.leaves,
					which(rows.to.update)
				)

		#remove first entry from all elements of list.of.node.names.to.prune[[]]
		list.of.node.names.to.prune[[list.index]] <- list.of.node.names.to.prune[[list.index]][-1]

		#create T/F vector of observations that change locations
		observations.that.move.T.F <- tree$where %in% which(row.names(tree$frame) %in% list.of.node.names.to.prune[[list.index]])
			tree$where[observations.that.move.T.F] <- which(rows.to.update)
			tree$y[observations.that.move.T.F] <- tree$frame$yval[rows.to.update]

		#update pruned row of tree$frame
			tree$frame$var[rows.to.update] <- "<leaf>"
			tree$frame$splits[rows.to.update,] <- matrix(c("",""),nrow=1)
	}

	#make a note of all row that will be removed
	rows.to.update <- row.names(tree$frame) %in% unique(unlist(list.of.node.names.to.prune))

	#go from old numbering system to new numbering system for 'tree$where'
	old.node.number <- which(!rows.to.update)
	if (length(old.node.number) > 1) {				
		for (new.node.number in 1:length(old.node.number)) {
			tree$where[tree$where == old.node.number[new.node.number]] <- new.node.number
		}
	} #where we have a stump, old and new numbering coincide

	#update the rest of tree
	tree$frame <- tree$frame[!rows.to.update,,drop=FALSE]		#remove pruned rows from tree$frame
	tree$details <- tree$details[!rows.to.update]			#remove pruned entires from tree$details

	#make sure leaves are indeed leaves, cant insert NULL into a list so have to create a new one
	all.leaves <- sapply(tree$details,is.null)
	all.leaves[new.leaves] <- TRUE
	new.details <- vector("list",length(all.leaves))
	new.details[!all.leaves] <- tree$details[!all.leaves]
	tree$details <- new.details

	#return updated tree object
	return(tree)
}

