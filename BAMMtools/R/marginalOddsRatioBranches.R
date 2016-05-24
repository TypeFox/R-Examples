
# returns a phylogenetic tree where 
# branch lengths are equal to the marginal 
# odds ratio (posterior : prior) for a given branch. 
# It is marginal in the sense that it is not independent 
# of values for other branches.

marginalOddsRatioBranches <- function(ephy, expectedNumberOfShifts) {
	
	tree_post <- marginalShiftProbsTree(ephy);
	tree_prior <- getBranchShiftPriors(as.phylo.bammdata(ephy), expectedNumberOfShifts);
	
	post_shift <- tree_post$edge.length;
	prior_shift <- tree_prior$edge.length;
 
	oddsratio <- post_shift / prior_shift;

	tree_post$edge.length <- oddsratio;
	return(tree_post);	
	
}



