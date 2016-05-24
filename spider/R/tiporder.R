tiporder <- function(phy, labels = TRUE){
	nn <- length(phy$tip.label) #How many tips on the tree?
	edge <- phy$edge
	nums <- rev(edge[edge[,2] %in% 1:nn, 2])
	if(labels == TRUE) phy$tip.label[nums] else nums
}
