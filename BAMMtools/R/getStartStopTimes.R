#############################################################
#
#	getStartStopTimes(....)
#
#	adds begin and end times (absolute time) to each edge of 
#	phylogenetic tree

getStartStopTimes <- function(phy){
 	if (is.ultrametric(phy)) {
	 	bmax <- max(branching.times(phy));
		bt <- bmax - branching.times(phy);
		begin <- bt[as.character(phy$edge[,1])];
		end <- begin + phy$edge.length;
		phy$begin <- as.numeric(begin);
		phy$end <- as.numeric(end);
		return(phy);
	}
	return( NU.branching.times(phy, "begin.end"));
}
