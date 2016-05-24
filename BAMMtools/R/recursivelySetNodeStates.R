#############################################################
#
#	recursivelySetNodeStates(....)
#
# Function to recursively assign states to nodes going from root-to-tip
#	takes args of tree plus vector of node event assignments
#	returns vector.
# phySV is tree with "statevec" component
# node is the focal node
# state is the state (probably "event" in this context)

recursivelySetNodeStates <- function(phySV, node, state) {
	
	phySV$statevec[phySV$edge[,2] == node] <- state;
	if (sum(phySV$edge[,1] == node) > 0){
		# node is internal
		dset <- phySV$edge[,2][phySV$edge[,1] == node];
		phySV <- recursivelySetNodeStates(phySV, dset[1], state);
		phySV <- recursivelySetNodeStates(phySV, dset[2], state);
	}
	
	return(phySV);
}
