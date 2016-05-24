#############################################################
#
#	getCladeRates(....)
#
#	mean clade-specific rates 
#		average of all branch-rates, but weighted by branch length
#		node.type: will compute rates only for clade descended from specified node with 'include'
#					will compute for all branches excluding a given clade, nodetype = 'exclude'
#		

getCladeRates <- function(ephy, node = NULL, nodetype='include', verbose=FALSE) {
	
	if (!'bammdata' %in% class(ephy)) {
		stop("Object ephy must be of class bammdata\n");
	}	
	
	if (is.null(node)) {
		nodeset <- ephy$edge[,2];
	} else if (!is.null(node[1]) & nodetype[1] == 'include' & length(node) == 1) {
		nodeset <- getDesc(ephy, node)$desc_set;
	} else if (!is.null(node[1]) & nodetype[1] == 'exclude' & length(node) == 1) {
		nodeset <- setdiff( ephy$edge[,2],  getDesc(ephy, node)$desc_set);
	} else if (!is.null(node[1]) & length(nodetype) == length(node) & length(node) > 1) {
		nodesets <- lapply(node, function(x) getDesc(ephy, x)$desc_set);
		Drop <- which(nodetype == 'exclude');
		nodeset_toRemove <- unique(unlist(lapply(nodesets[Drop], function(x) x)));
		Keep <- which(nodetype == 'include');
		nodeset_toKeep <- unique(unlist(nodesets[Keep]));
		nodeset <- setdiff(nodeset_toKeep, nodeset_toRemove);
		if (length(nodeset) == 0) {
			stop('Error: the combination of nodes and nodetypes has resulted in no remaining nodes!')
		}
	} else {
		stop('Error: Please make sure you have specified a nodetype for every node');
	}
	
	lambda_vector <- numeric(length(ephy$eventBranchSegs));
	mu_vector <- numeric(length(ephy$eventBranchSegs));
 	
 	weights <- 'branchlengths'
	
	for (i in 1:length(ephy$eventBranchSegs)) {
		if (verbose) {
			cat('Processing sample ', i, '\n');
		}
		esegs <- ephy$eventBranchSegs[[i]];
		esegs <- esegs[esegs[,1] %in% nodeset, ];
	
       	if (is.null(nrow(esegs))){
       		esegs <- t(as.matrix(esegs))
       	}	
		
		events <- ephy$eventData[[i]];
		events <- events[order(events$index), ];			
		
		# relative start time of each seg, to event:
		relsegmentstart <- esegs[,2] - ephy$eventData[[i]]$time[esegs[,4]];
		relsegmentend <- esegs[,3] - ephy$eventData[[i]]$time[esegs[,4]];
		lam1 <- ephy$eventData[[i]]$lam1[esegs[,4]];
		lam2 <-  ephy$eventData[[i]]$lam2[esegs[,4]];
		mu1 <-  ephy$eventData[[i]]$mu1[esegs[,4]];
		mu2 <-  ephy$eventData[[i]]$mu2[esegs[,4]];
 		
 		seglengths <- esegs[,3] - esegs[,2];	
		wts <- seglengths / sum(seglengths);
		lamseg <- timeIntegratedBranchRate(relsegmentstart, relsegmentend, lam1, lam2) / seglengths;
		museg <- timeIntegratedBranchRate(relsegmentstart, relsegmentend, mu1, mu2) / seglengths;
	
		lambda_vector[i] <- sum(lamseg * wts);
		mu_vector[i] <- sum(museg  * wts);
	}
	
	if (ephy$type == 'diversification') {
		return(list(lambda = lambda_vector, mu = mu_vector));
	}
	if (ephy$type == 'trait') {
		return(list(beta = lambda_vector));
	}
}
