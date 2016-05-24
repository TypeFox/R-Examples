recluster.check <-function(tree, tip) {
	ind <- match(tip, tree$tip.label)
	tree$tip.label <- tree$tip.label[ind]
	ind2 <- match(1:length(ind), tree$edge[, 2])
	tree$edge[ind2, 2] <- order(ind)
	tree
}
