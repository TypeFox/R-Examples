###############################################
#	Internal function called by plot.bammdata(...)
#
#
setPolarTreeCoords = function(phy,vtheta,rbf)
{
	phy <- getStartStopTimes(phy);
#	tH = max(branching.times(phy))
	tH <- max(phy$end);
	
	rootnd <- as.integer(phy$Nnode+2);
	ntip <- as.integer(phy$Nnode+1);
	anc <- as.integer(phy$edge[,1]);
	desc <- as.integer(phy$edge[,2]);
	nnode <- as.integer(dim(phy$edge)[1] + 1);
	
	vtheta <- as.numeric(vtheta*(pi/180));
	ths <- as.numeric((2*pi-vtheta)/phy$Nnode);	
	
	ndorder <- .C('postorder_tree_traverse', anc, desc, rootnd, nnode, integer(nnode), PACKAGE="BAMMtools")[[5]];
	ndorder <- as.integer(ndorder);
	L <- .C('setpolartreecoords', anc, desc, ndorder, ntip, rootnd, nnode, ths, numeric(nrow(phy$edge)*3), numeric(3), PACKAGE="BAMMtools");
	
	root <- matrix(L[[9]],nrow=1);
	theta <- matrix(L[[8]],nrow=nrow(phy$edge),ncol=3);
	theta <- rbind(root, theta);
	
	rb <- tH * rbf;
	x0 <- c(rb,rb+(phy$begin/tH))*cos(theta[,1]);
	y0 <- c(rb,rb+(phy$begin/tH))*sin(theta[,1]);
	x1 <- c(rb,rb+(phy$end/tH))*cos(theta[,1]);
	y1 <- c(rb,rb+(phy$end/tH))*sin(theta[,1]);
	ret <- cbind(x0,y0,x1,y1,theta[,1]);
	rownames(ret) <- c(phy$edge[1,1],phy$edge[,2]); colnames(ret) = c('x0','y0','x1','y1','theta');
	return(list(segs = ret, arcs = theta[,2:3]) );	
}
