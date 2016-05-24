#############################################################
#
#	getRecursiveSequence(....)
#
#	Private function, called by getEventDataDiversification

getRecursiveSequence = function(phy)
{
	rootnd = as.integer(phy$Nnode+2);
	anc = as.integer(phy$edge[,1]);
	desc = as.integer(phy$edge[,2]);
	ne = as.integer(dim(phy$edge)[1]);
	L = .C('setrecursivesequence', anc, desc, rootnd, ne, integer(ne+1),integer(ne+1));
	phy$downseq = as.integer(L[[5]]);
	phy$lastvisit = as.integer(L[[6]]);
	return(phy);
}