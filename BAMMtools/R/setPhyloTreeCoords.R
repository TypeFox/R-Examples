#######################################
#	Internal function called by plot.bammdata(...)
#
#
setPhyloTreeCoords = function(phy)
{
	phy <- getStartStopTimes(phy);
#	tH <- as.numeric(max(branching.times(phy)))
	tH <- max(phy$end);
		
	rootnd <- as.integer(phy$Nnode+2);
	ntip <- as.integer(phy$Nnode+1);
	anc <- as.integer(phy$edge[,1]);
	desc <- as.integer(phy$edge[,2]);
	nnode <- as.integer(dim(phy$edge)[1] + 1);
	bl <- as.numeric(phy$edge.length/tH);
	begin <- as.numeric(phy$end/tH);
	
	ndorder <- .C('postorder_tree_traverse', anc, desc, rootnd, nnode, integer(nnode), PACKAGE="BAMMtools")[[5]];
	ndorder <- as.integer(ndorder);
	
	L <- .C('setphylotreecoords', anc, desc, ndorder, begin, bl, ntip, rootnd, nnode, numeric(nrow(phy$edge)*4), numeric(nrow(phy$edge)*4), numeric(4), PACKAGE="BAMMtools");
	
	root <- matrix(L[[11]],nrow=1);
	xy <- matrix(L[[10]],nrow=nrow(phy$edge),ncol=4);
	bar <- matrix(L[[9]],nrow=nrow(phy$edge),ncol=4);
	
	xy <- rbind(c(xy[1,1],sum(xy[1:2,4])/2,xy[1,1],sum(xy[1:2,4])/2),xy);
	bar <- rbind(root,bar);
	rownames(xy) <- c(phy$edge[1,1],phy$edge[,2]); colnames(xy) = c('x0','y0','x1','y1');
	return(list (segs = xy, arcs = bar ) );	
}
