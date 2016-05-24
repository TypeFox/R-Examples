AC.node.trans <- function(new, old, index = FALSE){
	
	ntips <- length(new$tip.label)
	nodes_new <- (ntips + 1) : (ntips + new$Nnode)
	clades <- lapply(nodes_new, descendants, tree = new)
	name_clades <- function(x, phy){
		x <- phy$tip.label[x]
		x
	}
	clades <- lapply(clades, name_clades, phy = new)
	nodes_old <- noi(old, clades, monophyletic = TRUE)
	
	out <- cbind(nodes_new, nodes_old)
	if (index)
		out <- match(nodes_new, nodes_old)
	out
}


