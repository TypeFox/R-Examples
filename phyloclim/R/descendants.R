descendants <-
function(tree, node, internal = FALSE, string = FALSE){
	
	tips <- seq(along = tree$tip.label)
	x <- tree$edge[,2][tree$edge[,1] == node]
	repeat{
		xx <- x
		x <- sort(unique(c(x, tree$edge[,2][tree$edge[,1] %in% x])))
		if (identical(x, xx)) break
	}
	# return tip number if input is tip number:
	# -----------------------------------------
	if (length(x) == 0) x <- node
	if (!internal)
		x <- x[x %in% tips]
	if (string)
		x <- tree$tip.label[x]
	x
}