#############################################################
#
#	marginalShiftProbsTree(....)
#
#	Args: ephy	=	object of class 'bammdata'
#	
#	Returns:		a phylogenetic tree, but where each 
#	             	branch length (edge length) is equal to the
#					marginal probability of shift occuring 
#					on that particular branch. 
#							

marginalShiftProbsTree <- function(ephy) {
	
	if (!'bammdata' %in% class(ephy)) {
		stop("Object ephy must be of class bammdata\n");
	}
	
	shiftvec <- numeric(length(ephy$edge.length));
 	rootnode <- length(ephy$tip.label) + 1;
 
 
	for (i in 1:length(ephy$eventData)) {
		hasShift <- ephy$edge[,2] %in% ephy$eventData[[i]]$node;
		shiftvec[hasShift] <- shiftvec[hasShift] + rep(1, sum(hasShift));
	}
	
	shiftvec <- shiftvec / length(ephy$eventData);	
	
	newphy <- as.phylo.bammdata(ephy);
	newphy$edge.length <- shiftvec;
	return(newphy);
}
