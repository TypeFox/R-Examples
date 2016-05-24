HCLtree <- function(x){
	tmp.tree <- hclust(as.dist(1-x), "complete")
	return(tmp.tree)
}

