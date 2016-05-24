age.tips <-
function(tree){	
	times <- vector()
	if (class(tree)=="phylo"){
	tipdepth <- vector()
	for (j in (1:length(tree$tip))) {
		parent <- tree$edge[,1][tree$edge[,2]==j]
		tipdepthtemp <- branching.times.complete(tree)[as.character(parent)]
		tipdepthtemp <- c(j,round(tipdepthtemp - tree$edge.length[tree$edge[,2]==j],10))
		tipdepth <- rbind(tipdepth, tipdepthtemp)
		}
	times <- tipdepth
	}
	times
	}

