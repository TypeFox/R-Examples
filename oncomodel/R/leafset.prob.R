`leafset.prob` <-
function(leafset, y) {
	p <- y$p
        tree <- y$tree
        var.names <- y$var.names
  if (!is.tree(tree)) 
	stop("'tree' is not in suitable matrix format.")
  if (!is.vector(var.names)) 
      stop("'var.names' should be a vector.")
  if (!is.vector(leafset)) 
	stop("'leafset' is not in suitable vector format.")
  nodes <- sort(unique(as.vector(tree)))
  for (j in 1:length(nodes)) tree[which(tree == nodes[j])] <- j
  res <- as.list(c(var.names, rep("(", length(nodes) - length(var.names))))
  done <- rep(FALSE, length(nodes))
  done[1:length(var.names)] <- TRUE
  col.vec <- character(length(p))
  mat <- matrix(nrow=2,ncol=length(p))
  while (!all(done)) {
	waiting <- numeric()
      for (j in which(!done)) {
		if (all(done[tree[2, which(tree[1, ] == j)]])) 
            waiting <- c(waiting, j)
      }
      for (j in waiting) {
		edges <- which(tree[1, ] == j)
		for (k in 1:length(edges)) {
			col.vec[edges[k]] <- res[[tree[2, edges[k]]]]
			mat[1,edges[k]] <- edges[k]
			mat[2,edges[k]] <- as.numeric(format(p[edges[k]], digits=6))
			res[[j]] <- paste(res[[j]], res[[tree[2, edges[k]]]], sep = "")
			if (k == length(edges)) res[[j]] <- paste(res[[j]], ")", sep = "")
			else res[[j]] <- paste(res[[j]], ",", sep = "")
            }
            done[j] <- TRUE
        }
  }		# end of "while"-loop
  parent <- match(tree[1, ], tree[2, ])
  siblings <- list()
  for (j in 1:ncol(tree)) siblings[[j]] <- which(tree[1, ] == tree[1, j] & tree[2, ] != tree[2, j])
  siblings <- as.vector(siblings, mode="numeric")
  mat <- rbind(mat, (1-mat[2,]), parent, siblings)
  rownames(mat) <- c("column number", "Probability (respectively colNo.)", "'contra' probability",
			"parent (column number)", "sibling (column number)")
  colnames(mat) <- col.vec
  pos.leaves <- which(colnames(mat)%in%var.names)

  if (!all(is.element(which(colnames(mat)%in%leafset), pos.leaves)))
 	stop("There ist at least one invalid leafset value.")

#--- auxiliary functions
  #--- probability that non of the leaf events corresponding to the subtree occurs
  q_func <- function(node, pos.leaves=pos.leaves, mat=mat) {
	# node: root of the subtree
	# pos.leaves: position of the leaves in the data matrix
	# mat: data matrix
	if (!is.element(node, pos.leaves)) {
		children <- as.vector(which(mat[1,node]==mat[4,]))	# children: Nummer der Nachfolger
		if (is.element(children[1], pos.leaves) & is.element(children[2], pos.leaves))
			q_node <- mat[3,node]+mat[2,node]*mat[3,children[1]]*mat[3,children[2]]
		if (is.element(children[1], pos.leaves) & !is.element(children[2], pos.leaves))
			q_node <- mat[3,node]+mat[2,node]*mat[3,children[1]]*q_func(children[2], pos.leaves, mat)
		if (!is.element(children[1], pos.leaves) & is.element(children[2], pos.leaves))
			q_node <- mat[3,node]+mat[2,node]*mat[3,children[2]]*q_func(children[1], pos.leaves, mat)
		if (!is.element(children[1], pos.leaves) & !is.element(children[2], pos.leaves)) {
			q_children1 <- q_func(children[1], pos.leaves, mat)		# rekursiver Aufruf; alt: which(children[1]==mat[4,]), 
			q_children2 <- q_func(children[2], pos.leaves, mat)
			q_node <- mat[3,node]+mat[2,node]*q_children1*q_children2
		}
	} else q_node <- mat[3,node] # entspricht 1 - mat[2,node]
	return(q_node)
  }
  #--- calculation of part probability by bottom up procedure (path from leaf event to parent event)
  p_par <- function(node, mat=mat) {
	if (!is.na(mat[4, node])) {
		sibling <- mat[5, node]
		p_node <- mat[2,node]*q_func(sibling, pos.leaves, mat)*p_par(mat[4, node], mat)
	} else p_node <- mat[2,node]
	return(p_node)
  }
  #--- calculation of part probability by top down procedure (path from parent event to leaf event (child))
  p_child <- function(node.parent, k, nodes, mat=mat) {
	rest.nodes <- rev(unlist(nodes))
	p_vec <- numeric()
	if (k > 0) {
		sibling <- mat[5, rest.nodes[1]]	
		p_vec <- c(p_vec, q_func(sibling, pos.leaves, mat))
	}
	if (length(rest.nodes) <= 1) p_vec <- c(p_vec, mat[2, rest.nodes])
	else p_vec <- c(p_vec, mat[2, rest.nodes[1]], p_child(node.parent=rest.nodes[1], k=1, nodes=rev(rest.nodes[-1]), mat=mat))
	return(p_vec)
  }
  #--- probability for two branches
  p_2_branches <- function(list, gem.par, pos, mat=mat) {
	# list: path to leaf event
	# gem.par: joint parent event
	p_vec <- numeric()
	for (j in 1:2) {	# length(leafset)) {
		rest.vec <- vector()
		k <- 0
		rest.vec <- list[[j]][-which(list[[j]] %in% gem.par)]
		if (length(rest.vec) <= 2)
			p_vec <- c(p_vec, p_child(node.parent=min(gem.par), k=k, nodes=rest.vec, mat=mat))
		else {
			p_vec <- c(p_vec, p_child(node.parent=min(gem.par), k=k, nodes=rest.vec, mat=mat))
		}	#mat[2, rev(rest.vec)[1]], 
	}
	return(list(rest.vec, p_vec))
  }

#--- calculation of probability
  n <- length(unique(leafset))
  leafset <- unique(leafset)
  pos <- match(leafset, colnames(mat))
  p_vec <- numeric()
  if (n <= 1) {
  # one change/leave
	sibling <- mat[5, pos]
	p_vec <- c(mat[2, pos], q_func(sibling, pos.leaves, mat), p_par(mat[4, pos], mat))
  } else {
  # two or more changes
	list <- list()
	for (i in 1:length(leafset)) {
		index <- paste(i, ": path to ", leafset[i], sep="")
		list[[index]] <- grep(leafset[i], colnames(mat), fixed=TRUE)
              }
	# all leafsets => from/ex root
	if (n >= length(pos.leaves)) p_vec <- mat[2,]
	else {
	  repeat {
		nr1 <- nr2 <- numeric()
		for (i in 1:(length(pos)-1)) {
			for (j in (i+1):length(pos)) {
				nr1 <- c(nr1, i)
				nr2 <- c(nr2, j)
		}	}
		gem.par.list <- list()
		index_vec <- vector()
		length_vec <- vector()
		for (k in 1:length(nr1)) {
			index <- paste(nr1[k], nr2[k], sep=":")
			index_vec <- c(index_vec, index)
			gem.par.list[[index]] <- intersect(list[[nr1[k]]], list[[nr2[k]]])
			length_vec <- c(length_vec, length(gem.par.list[[k]]))
		}
		index.pos <- which(max(length_vec)==length_vec)
		list.nr <- as.numeric(unlist(strsplit(index_vec[index.pos],":")))
		# siblings
		if (mat[5,list[[list.nr[1]]][1]] == mat[1,list[[list.nr[2]]][1]] & mat[5,list[[list.nr[2]]][1]] == mat[1,list[[list.nr[1]]][1]]) {
			p_vec <- c(p_vec, mat[2,pos[list.nr[1]]], mat[2,pos[list.nr[2]]])
			list[[list.nr[1]]] <- list[[list.nr[1]]][-1]
			pos[list.nr[1]] <- list[[list.nr[1]]][1]
			list <- list[-list.nr[2]]
			pos <- pos[-list.nr[2]]
		} else {
			sub.gem.par <- intersect(list[[list.nr[1]]],list[[list.nr[2]]])
			sub.list <- list(list[[list.nr[1]]],list[[list.nr[2]]])
			res <- p_2_branches(list=sub.list, gem.par=sub.gem.par, pos=pos[list.nr], mat=mat)
			p_vec <- c(p_vec, res[[2]])
			# delete all nodes up to first shared node
			list[[list.nr[1]]] <- list[[list.nr[1]]][-(1:length(setdiff(list[[list.nr[1]]],sub.gem.par)))]
			pos[list.nr[1]] <- min(sub.gem.par)
			list <- list[-list.nr[2]]
			pos <- pos[-list.nr[2]]
		}
		if (length(pos) < 2) {
			if (pos != ncol(mat)) p_vec <- c(p_vec, p_par(pos, mat))
			break
		}
	  }	# end of "repeat"-loop
	}
  }		# end of "if (n > 1)" (else)
  p_return <- prod(p_vec)
  return(p_return)
}

