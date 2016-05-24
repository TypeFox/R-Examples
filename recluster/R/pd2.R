pd2 <- function (samp, tree, include.root = TRUE) {
	if (is.null(tree$edge.length)) {
		stop("Tree has no branch lengths, cannot compute pd")
		}
	species <- colnames(samp)
	tree	<- node.age(tree)
	PDout <- apply(samp,1, function(x) {
		present <- species[x > 0]
		treeabsent <- tree$tip.label[which(!(tree$tip.label %in%present))]
		if (length(present) == 0) {
			PD <- 0
		}
		else if (length(present) == 1) {
			if (!is.rooted(tree) || !include.root) {
				warning("Rooted tree and include.root=TRUE argument required to calculate PD of single-species sampunities. Single species sampunity assigned PD value of NA.")
				PD <- NA
			}

			else {
				PD <- tree$ages[which(tree$edge[, 2] ==
				which(tree$tip.label == present))]
			}
		}	
		else if (length(present[!duplicated(present)]) == 1) {
			if (!is.rooted(tree) || !include.root) {
				warning("Rooted tree and include.root=TRUE argument required to calculate PD of single-species sampunities. Single species sampunity assigned PD value of NA.")
				PD <- NA
			}
			else {
				PD <- tree$ages[which(tree$edge[, 2] ==
				which(tree$tip.label == present[1]))]
			}
		}	
		else if (length(treeabsent) == 0) {
			PD <- sum(tree$edge.length)
		}
		else {
			sub.tree <- drop.tip(tree, treeabsent)
			if (include.root) {
				if (!is.rooted(tree)) {
					stop("Rooted tree required to calculate PD with include.root=TRUE argument")
				}
				sub.tree <- node.age(sub.tree)
				sub.tree.depth <- max(sub.tree$ages)
				orig.tree.depth <- max(tree$ages[which(tree$edge[,2] %in% which(tree$tip.label %in% present))])
				PD <- sum(sub.tree$edge.length) + (orig.tree.depth - sub.tree.depth)
			}
			else {
				PD <- sum(sub.tree$edge.length)
			}
		}	 
		SR <- length(present)	
		PDout <- c(PD,SR)
	} )				 
	PDout <- t(PDout)
	rownames(PDout) <- rownames(samp)
	colnames(PDout) <- c("PD","SR")	 
	return(PDout)	
}
