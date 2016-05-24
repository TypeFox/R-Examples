addBAMMshifts = function(ephy, index = 1, method = 'phylogram', cex=1, pch=21, col=1, bg=2, msp = NULL, shiftnodes = NULL, par.reset=TRUE) {
	if (!'bammdata' %in% class(ephy)) stop("Object ephy must be of class bammdata");
	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv);
	
	if (par.reset){
		op <- par(no.readonly = TRUE);
		par(lastPP$pp);		
	}

	if (length(ephy$eventData) == 1){
		index <- 1;
	}
	
	if (is.null(shiftnodes))
		shiftnodes <- getShiftNodesFromIndex(ephy, index)
	isShift <- ephy$eventData[[index]]$node %in% shiftnodes;
	times <- ephy$eventData[[index]]$time[isShift];	
	if (!is.null(msp)) {
		cex <- 0.75 + 5 * msp$edge.length[msp$edge[,2] %in% shiftnodes];
	}
	
	if (method == 'phylogram') {
		###  obsolete b/c plot.bammdata no longer scales each axis to a max of 1. now behaves like plot.phylo
		# if (max(lastPP$xx) <= 1) {
		# 	XX <- times/max(branching.times(as.phylo.bammdata(ephy)));
		# } else {
		# 	XX <- times;
		# }
		XX <- times;
		YY <- lastPP$yy[shiftnodes];
	} else if (method == 'polar') {
		rb <- lastPP$rb;
		XX <- (rb+times/max(branching.times(as.phylo.bammdata(ephy)))) * cos(lastPP$theta[shiftnodes]);
		YY <- (rb+times/max(branching.times(as.phylo.bammdata(ephy)))) * sin(lastPP$theta[shiftnodes]);		
	}	
	points(XX,YY,pch=pch,cex=cex,col=col,bg=bg);
	if (par.reset) {
		par(op);		
	}
}