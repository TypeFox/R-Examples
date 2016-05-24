

testTimeVariableBranches <- function(ephy, prior_tv = 0.5, return.type = 'posterior'){

	# return.type = bayesfactor or posterior

	TOL <- .Machine$double.eps * 10;	
	
	esum <- numeric(nrow(ephy$edge));
	
	em <- ephy$edge;
	
	for (i in 1:length(ephy$eventBranchSegs)){
		
		segmat <- ephy$eventBranchSegs[[i]];
		rv <- ephy$eventData[[i]]$lam2[segmat[,4]];
		rv <- as.numeric(abs(rv) > TOL);
		
		for (k in 1:nrow(ephy$edge)){
			
			esum[k] <- esum[k] + rv[which(segmat[,1] == em[k,2])[1]];
			
			
		}
	}
	prob.rv <- esum / length(ephy$eventData);
	
	bl <- prob.rv;
	if (return.type == "bayesfactor"){
		prob.null <- 1 - prob.rv;
		bl <- (prob.rv / prob.null) * ((1 - prior_tv) / prior_tv);
		
	}else if (return.type != "posterior"){
		stop("Invalid return type\n");
		
	}
	
	
	newphy <- list(edge=ephy$edge, edge.length = bl, Nnode = ephy$Nnode, tip.label=ephy$tip.label);
	
	class(newphy) <- 'phylo';
	newphy;
}



