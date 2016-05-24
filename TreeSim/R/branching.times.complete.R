branching.times.complete <-
function(tree){
    if (class(tree) != "phylo") 
        stop("object \"tree\" is not of class \"phylo\"")
    n <- length(tree$tip.label)
    N <- dim(tree$edge)[1]
    xx <- numeric(tree$Nnode)
    interns <- which(tree$edge[, 2] > n)
    for (i in interns) xx[tree$edge[i, 2] - n] <- xx[tree$edge[i, 
        1] - n] + tree$edge.length[i]
    depthtemp <- xx[tree$edge[, 1] - n] + tree$edge.length
    depth <- max(depthtemp)
    xx <- depth - xx
    names(xx) <- if (is.null(tree$node.label)) 
        (n + 1):(n + tree$Nnode)
    else tree$node.label
      xx
	}

