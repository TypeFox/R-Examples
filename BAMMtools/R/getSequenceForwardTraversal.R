#############################################################
#
#	getSequenceForwardTraversal(....)
#
#	Private function, called by getRecursiveSequence


getSequenceForwardTraversal <- function(phy, node){
	
	if (node <= length(phy$tip.label)){
		#phy$downseq <- c(phy$downseq, node);
	}else{
		dset <- phy$edge[,2][phy$edge[,1] == node];
		phy$downseq <- c(phy$downseq, dset[1]);
		phy <- getSequenceForwardTraversal(phy, dset[1]);
		phy$downseq <- c(phy$downseq, dset[2]);
		phy <- getSequenceForwardTraversal(phy, dset[2]);	
	}
	
	phy$lastvisit[node] <- phy$downseq[length(phy$downseq)];
 
	return(phy);
}
