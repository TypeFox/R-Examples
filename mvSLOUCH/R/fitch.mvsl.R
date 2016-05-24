## Function from the slouch package. Author : Jason Pienaar
'fitch.mvsl'<-function (phyltree, niche, deltran = FALSE, acctran = FALSE, 
root = NULL) 
{
## ----- Added Krzysztof Bartoszek ----------------------------------------
    tree.data<-.ouch2slouch.mvsl(phyltree)
    niche<-c(rep(NA,phyltree@nnodes-phyltree@nterm),as.character(niche))
## ------------------------------------------------------------------------
    
    niche <- as.factor(niche)
    anc.char <- as.character(tree.data$ancestor)
    anc.num <- as.numeric(anc.char)
    anc.num[1] <- 0
    node.char <- as.character(tree.data$nodes)
    node.num <- as.numeric(node.char)
    niche.num <- as.numeric(niche)
    niches <- levels(niche)
    num.code <- 1:length(levels(niche))
    translator <- cbind(niches, num.code)
    obs.char <- unique(niche.num[!is.na(niche.num)])
    N <- length(node.num)
    tree.matrix <- matrix(data = c(anc.num, node.num, niche.num), 
						  ncol = 3, )
    colnames(tree.matrix) <- c("Ancestors", "Nodes", "States")
    Nint <- max(tree.matrix[, 1])
    Nnodes <- max(tree.matrix[, 2])
    downpass.node.states <- list()
    cost <- 0
    for (i in 1:Nnodes) {
        downpass.node.states[[i]] <- tree.matrix[i, 3]
    }
    check.node.order <- tree.matrix[, 2] - tree.matrix[, 1]
    if (any(check.node.order < 0)) 
	traversal <- 2:Nint
    else traversal <- Nint:1
    for (i in traversal) {
        children <- which(tree.matrix[i, 2] == tree.matrix[, 
						  1])
        if (length(children) == 2) {
            child1.state <- downpass.node.states[[children[1]]]
            child2.state <- downpass.node.states[[children[2]]]
            downpass.node.states[[i]] <- intersect(child1.state, 
												   child2.state)
        }
        if (length(children) > 2) {
            child1.state <- downpass.node.states[[children[1]]]
            child2.state <- downpass.node.states[[children[2]]]
            child3.state <- downpass.node.states[[children[3]]]
            tmp <- intersect(child1.state, child2.state)
            downpass.node.states[[i]] <- intersect(tmp, child3.state)
        }
        if (length(downpass.node.states[[i]]) == 0) {
            if (length(children) == 2) {
                downpass.node.states[[i]] <- union(child1.state, 
												   child2.state)
                cost = cost + 1
            }
            if (length(children) > 2) {
                tmp <- union(child1.state, child2.state)
                downpass.node.states[[i]] <- union(tmp, child3.state)
                cost = cost + 1
            }
        }
    }
    if (any(check.node.order < 0)) {
        children <- which(tree.matrix[1, 2] == tree.matrix[, 
						  1])
        if (length(children) == 2) {
            child1.state <- downpass.node.states[[children[1]]]
            child2.state <- downpass.node.states[[children[2]]]
            downpass.node.states[[1]] <- intersect(child1.state, 
												   child2.state)
        }
        if (length(children) > 2) {
            child1.state <- downpass.node.states[[children[1]]]
            child2.state <- downpass.node.states[[children[2]]]
            child3.state <- downpass.node.states[[children[3]]]
            downpass.node.states[[1]] <- intersect(c(child1.state, 
													 child2.state), child3.state)
        }
        if (length(downpass.node.states[[1]]) == 0) {
            if (length(children) == 2) {
                downpass.node.states[[1]] <- union(child1.state, 
												   child2.state)
                cost = cost + 1
            }
            if (length(children) > 2) {
                tmp <- union(child1.state, child2.state) ## changed from   tmp <- union(child1.state, child.state) by Krzysztof Bartoszek
                downpass.node.states[[1]] <- union(tmp, child3.state)
                cost = cost + 1
            }
        }
    }
    if (any(check.node.order < 0)) 
	traversal <- 2:Nint
    else traversal <- Nint:2
    pre.traversal <- rev(traversal)
    finalpass.node.states <- downpass.node.states
    if (!is.null(root)) 
	finalpass.node.states[[1]] <- root
    if (length(finalpass.node.states[[1]]) >= 2) {
        message("There is an ambiguity at the root node, as given below")
        print(translator[as.numeric(finalpass.node.states[[1]])])
        message("Set this to one of the ambiguous states using root = state (state needs to be in quotation marks) in the function call before attempting the deltran or acctran reconstructions")
    }
    for (i in pre.traversal) {
        parent <- tree.matrix[i, 2]
        ancestor <- tree.matrix[i, 1]
        children <- which(tree.matrix[, 1] == tree.matrix[i, 
						  2])
        finalpass.node.states[[parent]] <- intersect(downpass.node.states[[parent]], 
													 finalpass.node.states[[ancestor]])
        if (setequal(finalpass.node.states[[parent]], finalpass.node.states[[ancestor]]) == 
            FALSE) {
            if (length(children) == 2) {
                child1.state <- downpass.node.states[[children[1]]]
                child2.state <- downpass.node.states[[children[2]]]
                if (length(intersect(child1.state, child2.state)) != 
					0) {
					finalpass.node.states[[parent]] <- union(downpass.node.states[[parent]], 
															 (intersect(finalpass.node.states[[ancestor]], 
																		union(child1.state, child2.state))))
                }
                if (length(intersect(child1.state, child2.state)) == 
					0) {
					finalpass.node.states[[i]] <- union(downpass.node.states[[parent]], 
														finalpass.node.states[[ancestor]])
                }
                if (deltran == TRUE) {
					if (length(finalpass.node.states[[i]]) >= 2) 
                    finalpass.node.states[[i]] = finalpass.node.states[[ancestor]]
                }
                if (acctran == TRUE) {
					if (length(finalpass.node.states[[i]]) >= 2) 
                    finalpass.node.states[[i]] = setdiff(finalpass.node.states[[parent]], 
														 finalpass.node.states[[ancestor]])
                }
            }
            if (length(children) > 2) {
                child1.state <- downpass.node.states[[children[1]]]
                child2.state <- downpass.node.states[[children[2]]]
                child3.state <- downpass.node.states[[children[3]]]
                if (length(intersect(child3.state, intersect(child1.state, 
															 child2.state))) != 0) {
					finalpass.node.states[[parent]] <- union(downpass.node.states[[parent]], 
															 (intersect(finalpass.node.states[[ancestor]], 
																		union(union(child3.state, child1.state), 
																			  child2.state))))
                }
                if (length(intersect(child3.state, intersect(child1.state, 
															 child2.state))) == 0) {
					finalpass.node.states[[i]] <- union(downpass.node.states[[parent]], 
														finalpass.node.states[[ancestor]])
                }
            }
        }
    }
    dat <- list(treematrix = tree.matrix, Final.states = finalpass.node.states, 
				downpass.cost = cost, niche.code = translator)
    N <- length(dat$Final.states)
    count <- 0
    ambig <- NA
    for (i in 1:N) {
        if (length(dat$Final.states[[i]]) >= 2) {
            count = count + 1
            ambig <- c(ambig, i)
        }
    }
    if (count != 0) {
        message("The number of ambiguous nodes are:")
        print(count)
        message(" and they are at nodes:")
        print(ambig[-1])
        message("Try deltran OR acctran = TRUE in the function call in order to implement a delayed or accelerated character transformation")
    }
    f.states <- .make.states.mvsl(dat)
    return(f.states)
}
