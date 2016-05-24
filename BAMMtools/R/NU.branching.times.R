
# Recursively compute branching times for phylogenetic tree.
#	Allows for non-ultrametric (fossil) trees.
#	Two return types:
#		return.type == 'bt'
#			returns traditional ape-like branching times (branching.times)
#		return.type == 'begin.end'
#			returns phylogeny with begin and end vectors
#				exactly like getStartStopTimes
#
#		This is slower than ape::branching.times, on account of
#			the recursion. Could recode in C for speed.
#
NU.branching.times <- function(phy, return.type = 'bt'){
	
	if (!is.binary.tree(phy)){
		stop("error. Need fully bifurcating (resolved) tree\n")
	}
	
	phy$begin <- rep(0, nrow(phy$edge)) 
	phy$end <-  rep(0, nrow(phy$edge))
	
	# Do it recursively
	
	fx <- function(phy, node){
	
		cur.time <- 0
		root <- length(phy$tip.label) + 1
		if (node > root){
			cur.time <- phy$end[which(phy$edge[,2] == node)]	
		}
	
		dset <- phy$edge[,2][phy$edge[,1] == node]
	
		i1 <- which(phy$edge[,2] == dset[1])
		i2 <- which(phy$edge[,2] == dset[2])
	
		phy$end[i1] <- cur.time + phy$edge.length[i1]
		phy$end[i2] <- cur.time + phy$edge.length[i2]

		if (dset[1] > length(phy$tip.label)){
			phy$begin[phy$edge[,1] == dset[1]] <- phy$end[i1]
			phy <- fx(phy, node = dset[1])	
			
		}
		if (dset[2] > length(phy$tip.label)){
			phy$begin[phy$edge[,1] == dset[2]] <- phy$end[i2]
			phy <- fx(phy, node = dset[2])
			
		}
		
		return(phy)
		
	}
	
	phy <- fx(phy, node = length(phy$tip.label) + 1)
	
	if (return.type == 'bt'){
		
		maxbt <- max(phy$end)
		nodes <- (length(phy$tip.label) + 1):(2*length(phy$tip.label) - 1)
		bt <- numeric(length(nodes))
		names(bt) <- nodes
		for (i in 1:length(bt)){
			tt <- phy$begin[phy$edge[,1] == nodes[i]][1]
			bt[i] <- maxbt - tt
		}
		
		return(bt)
		
	}else if (return.type == 'begin.end'){
		return(phy)		
	}
	
}


