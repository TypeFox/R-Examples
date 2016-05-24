


# Drop all nodes from event data with marginal probs < 0.05
# test is unique
#  if so, add to list
# Returns:
#	$marg.probs = marginal probs for nodes
#	$marginal_odds_ratio = branch-specific (marginal) posterior:prior odds ratios associated with 1 or more shifts
#	$shifts = unique shift sets
#	$samplesets = list of sample indices that reduce to each of the unique shift sets
#	$frequency = vector of frequencies of each shift configuration
#	$threshold =  (marginal) posterior:prior odds ratio threshold for shifts
#	
#	Results are sorted by frequency. 
#	$frequency[1] gives the most common shift config sampled
#	$shifts[[1]] gives the corresponding node indices for that configuration
#	$samplesets[[1]] gives the indices of samples with this configuration



distinctShiftConfigurations <- function(ephy, expectedNumberOfShifts, threshold, ... ) {
  
	or <- marginalOddsRatioBranches(ephy, expectedNumberOfShifts)
	
	mm <- marginalShiftProbsTree(ephy);

	goodnodes <- or$edge[,2][or$edge.length >= threshold];

	xlist <- list();
	for (i in 1:length(ephy$eventData)) {
		xlist[[i]] <- intersect(goodnodes, ephy$eventData[[i]]$node);
	}

	ulist <- list();
	treesets <- list();
	
	ulist[[1]] <- xlist[[1]];
	treesets[[1]] <- 1;
	
	for (i in 2:length(xlist)) {
		lx <- length(ulist);
		#cat(lx, '\n')
		for (k in 1:lx) {
			if (areShiftSetsEqual(ulist[[k]], xlist[[i]])){
				treesets[[k]] <- c(treesets[[k]], i);
				break;	
			} else {
				if (k == length(ulist)){
					xlen <- length(ulist);
					ulist[[xlen + 1]] <- xlist[[i]];
					treesets[[xlen + 1]] <- i;
				}
			}
		}
	}
	
	freqs <- unlist(lapply(treesets, length));
	freqs <- freqs / sum(freqs);
	
	ord <- order(freqs, decreasing=TRUE);
	
	obj <- list();
	obj$marg.probs <- mm$edge.length;  
	names(obj$marg.probs) <- mm$edge[,2]; 
	obj$marginal_odds_ratio <- or$edge.length;
	names(obj$marginal_odds_ratio) <- or$edge[,2];
	obj$shifts <- ulist[ord]; 
	obj$samplesets <- treesets[ord];
	obj$frequency <- freqs[ord];
	obj$coreshifts <- goodnodes;
	obj$threshold <- threshold;

	class(obj) <- 'bammshifts';
	
	return(obj);
}






