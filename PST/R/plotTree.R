## ==============
## Plotting nodes
## ==============
plotTree <- function(x1, x2, subtree, seglist, nPar, ePar, horiz = FALSE, gratio, max.level, cex, nc, cpal) {

	## Retrieving requested attributes
	cpal <- Xtract("cpal", nPar, default = NULL)
	depth <- subtree@order

	children <- which.child(subtree)

	## inner <- !length(children)==0 && x1 != x2 && !depth==max.level
	inner <- !length(children)==0 && !depth==max.level
	yTop <- subtree@order

	bx <- plotNodeLimit(x1, x2, subtree, max.level)
	xTop <- bx$x

	## If there are children, plotting edges and child nodes
	if (inner) {
		plotEdge(x1, x2, subtree, ePar, horiz, gratio, max.level, cex, nc, cpal)

		## Plotting edges and child nodes
		## Selecting non null child nodes only
		for (k in children) {
			child <- subtree[[k]]
			idx <- which(k==children)

			## Plotting the subtree
			plotTree(bx$limit[idx], bx$limit[idx+1], subtree = child, seglist=seglist, nPar=nPar, ePar=ePar, 
                		horiz=horiz, gratio=gratio, max.level=max.level, cex=cex, nc=nc)
		}
	} 

	plotNode(x1, x2, subtree, seglist, nPar, horiz, gratio, max.level, cex, nc, cpal)

}

