#############################################################
#
#	cumulativeShiftProbsTree(....)
#
#	Args: ephy	=	object of class 'bammdata'
#	
#	Returns:		a phylogenetic tree, but where each 
#	             	branch length (edge length) is equal to the
#					cumulative probability of a shift somewhere 
#					between the focal branch and the root of the 
#					tree. The branch length itself does not tell 
#					you where the shifts occur, but they tell 
#					you which clades/lineages have diversification
#					dynamics that are decoupled from the root of the tree 
#							

cumulativeShiftProbsTree <- function(ephy) {
	
	if (!'bammdata' %in% class(ephy)) {
		stop("Object ephy must be of class bammdata\n");
	}	
			
	shiftvec <- numeric(length(ephy$edge.length));
 	rootnode <- length(ephy$tip.label) + 1;
 
	for (i in 1:length(ephy$eventData)) {
		snodes <- unique(ephy$eventBranchSegs[[i]][,1][ephy$eventBranchSegs[[i]][,4] != 1]);
		hasShift <- ephy$edge[,2] %in% snodes;
		shiftvec[hasShift] <- shiftvec[hasShift] + rep(1, sum(hasShift));
	}	
	
	shiftvec <- shiftvec / length(ephy$eventData);	
	newphy <- as.phylo.bammdata(ephy);
	newphy$edge.length <- shiftvec;
	return(newphy);
}
