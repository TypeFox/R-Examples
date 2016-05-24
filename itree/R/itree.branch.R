#ALG: internal function that itree uses with slight modifications
# from rpart to deal with highlighting.

#SCCS @(#)rpart.branch.s	1.2 01/25/97
#
# Compute the "branches" to be drawn for an rpart object
#
itree.branch <- function(x, y, node, branch) {
	#ALG 9/7/2012: x and y are the coordinates of the nodes
	#in order of the tree's 'frame' object.

    if (missing(branch)) {
    	#6/24. made to match current rpart version in which envir!=Global environment
		if (exists(parms <-paste(".itree.parms", dev.cur(), sep="." ),
                   envir=itree_env)) {
#	    	parms <- get(parms, frame=0)
            parms <- get(parms, envir=itree_env) #ALG: same here, change to itree_env
            branch <- parms$branch
	    }
			else branch <- 0
        }

    # Draw a series of horseshoes, left son, up, over, down to right son
    #   NA's in the vector cause lines() to "lift the pen"
    is.left <- (node%%2 ==0)        #left hand sons
    node.left <- node[is.left]
    parent <- match(node.left/2, node)
    sibling <- match(node.left+1, node)
    temp <- (x[sibling] - x[is.left])*(1-branch)/2

    xx <- rbind(x[is.left], x[is.left]+ temp,
                x[sibling]- temp, x[sibling], NA)
    yy <- rbind(y[is.left], y[parent], y[parent], y[sibling], NA)

    #ALG 9/7/2012: xx and yy are matrices where each column represents
    #a left-node. the ith column gives the coordinates for drawing the horseshoe
    #starting at the ith left node.
    # In addition to xx and yy, I return a list of the left-node, right-node, and parent node
    # rownumbers in 'frame' corresponding to each column of xx & yy.
    list(x=xx, y=yy, leftnode.row=match(node.left,node), rightnode.row=sibling, parentnode.row=parent)
    }
