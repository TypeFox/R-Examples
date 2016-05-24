### subsetEventData

subsetEventData <- function(ephy, index) {
	
	if (class(ephy) != 'bammdata') {
		stop("Object ephy must be of class bammdata\n");
	}
	
	nsamples <- length(ephy$eventData);

	ss <- which((1:nsamples) %in% index);
	badsubset <- setdiff(index, 1:nsamples);
	
	if (length(badsubset) > 0) {
		cat("Bad call to BAMMtools::subsetEventData\n");
		cat("You have << ", nsamples, " >> samples in your data object \n");
		stop("Attempt to access invalid samples. Check index.")
	}
	
	obj <- list();
	obj$edge <- ephy$edge;
	obj$Nnode <- ephy$Nnode;
	obj$tip.label <- ephy$tip.label;
	obj$edge.length <- ephy$edge.length;
	obj$begin <- ephy$begin;
	obj$end <- ephy$end;
	obj$downseq <- ephy$downseq;
	obj$lastvisit <- ephy$lastvisit;
	
	obj$numberEvents <- ephy$numberEvents[ss];
	obj$eventData <- ephy$eventData[ss];	
	obj$eventVectors <- ephy$eventVectors[ss];
	obj$tipStates <- ephy$tipStates[ss]
	obj$tipLambda <- ephy$tipLambda[ss];
	obj$tipMu <- ephy$tipMu[ss];
	obj$eventBranchSegs <- ephy$eventBranchSegs[ss];
	
	obj$meanTipLambda <- ephy$meanTipLambda;
	obj$meanTipMu <- ephy$meanTipMu;
	
	obj$type <- ephy$type;
	
	class(obj) <- 'bammdata';
	attributes(obj)$order = attributes(ephy)$order;	
	return(obj);
}
