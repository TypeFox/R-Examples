#############################################################
#
#	getShiftNodesFromIndex (....)
#
#	Args: ephy	=	object of class 'bammdata'
#		  index =   the index of the sample you wish to view, e.g.,
#					 if index = 5, this will give you the nodes subtending
#                    all branches with rate shifts for the 5th sample
#					 from the posterior in your 'bammdata' object.								 
#
#	
#	Returns: 		- a vector of the nodes where shifts occurred, excluding the root.
#					Note: if NO shifts occured, this will return a 
#							numeric object of length zero
# 

getShiftNodesFromIndex <- function(ephy, index) {
	
	if (index > length(ephy$eventData)) {
		cat("Error. Attempting to access non-existent element from 'bammdata' object\n");
		cat("You have << ", length(ephy$eventData), " >>> samples in your 'bammdata' object\n");
		cat("Value of index must be no greater than the number of samples\n\n");
		stop();
	}
	root <- length(ephy$tip.label) + 1;
	nodes <- ephy$eventData[[index]]$node;
	nodes <- nodes[nodes != root];

	return(nodes);
}
