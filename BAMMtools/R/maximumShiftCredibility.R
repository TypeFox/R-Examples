#############################################################
#
#	maximumShiftCredibility(....)
#
#	Args: ephy	=	object of class 'bammdata'
#	
#
#	maximize		=	'sum', = sum of branch probabilities for each tree
#						'product', = product of branch probabilities for each tree
#	
#	Returns: 		- bestconfigs: a list of length equal the number of 
#						unique shift configurations in the maximum shift
#						credibility set. Each element is a vector of sample
#						indices from the 'bammdata' object with identical
#						shift configurations.
#					  
#					- A vector of optimality scores for all other samples 
#						in posterior from the 'bammdata' object.
#
#					- sampleindex: a representative index for samples from each 
#						set of unique shift configurations. The length of this vector
#						is equal to the length of the bestconfigs list. If this vector was
#						sampleindex = c(2, 20, 50), this would mean that there are 3 distinct
#						sets of shift configurations with equal credibility under the optimality 
#						criterion. More commonly, a single shift configuration will be dominant, and 
#						although the length of bestconfigs[[1]] may be greater than 1, the sampleindex
#						vector will contain a single representative event from that set.
#
#	See example file.
#   This is analogous to the maximum clade credibility tree from a 
#		Bayesian phylogenetic analysis.

maximumShiftCredibility <- function(ephy, maximize = 'product') {

	if (!'bammdata' %in% class(ephy)) {
		stop("Object ephy must be of class bammdata\n");
	}			
	
	probvec <- numeric(length(ephy$eventData));
	
	mtree <- marginalShiftProbsTree(ephy);
	
	#mtree$edge.length[mtree$edge.length < threshold] <- 0;
	
	px <- mtree$edge.length;
	
	ttx <- table(ephy$numberEvents) / length(ephy$numberEvents);
	
	
	for (i in 1:length(ephy$eventData)) {
		
		# posterior probabilities here:
		proc_prob <- ttx[as.character(ephy$numberEvents[i])]; 
		
		hasShift <- ephy$edge[,2] %in% ephy$eventData[[i]]$node;
		branchprobs <- (hasShift)*px  + (!hasShift)*(1 - px) ;
		if (maximize == 'product') {
			probvec[i] <- log(proc_prob) + sum(log(branchprobs));
		} else if (maximize == 'sum') {
			probvec[i] <- proc_prob * sum(branchprobs);
		} else {
			stop("Unsupported optimize criterion in maximumShiftCredibilityTree");
		}
	}
	
	best <- which(probvec == max(probvec));
	
	# Now test for multiple trees with same log-prob:
	bestconfigs <- list();
		
	index <- 0;	
	while (length(best) > 0) {
		index <- index + 1;	
		lv <- logical(length = length(best));
		for (i in 1:length(best)) {
			lv[i] <- areEventConfigurationsIdentical(ephy, best[1], best[i]);
		}
		bestconfigs[[index]] <- best[lv];
		best <- best[!lv];
	}
	
	sampleindex <- numeric(length(bestconfigs));
	for (i in 1:length(bestconfigs)) {
		sampleindex[i] <- bestconfigs[[i]][1];
	}
	
	obj <- list();
	obj$bestconfigs <- bestconfigs;
	obj$scores <- probvec;
	obj$optimalityType = maximize;
 	obj$sampleindex <- sampleindex;
	return(obj);
}
