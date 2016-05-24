
# Feb 28 2014
credibleShiftSet <- function(ephy, expectedNumberOfShifts, threshold = 5, set.limit = 0.95, ...){
	
	prior <- getBranchShiftPriors(ephy, expectedNumberOfShifts)
	
	dsc <- distinctShiftConfigurations(ephy, expectedNumberOfShifts, threshold);
	cfreq <- cumsum(dsc$frequency);
	cut <- min(which(cfreq >= set.limit));
	nodeset <- NULL;
 	
 	shiftnodes <- dsc$shifts[1:cut];
	indices <- dsc$samplesets[1:cut];
	frequency <- dsc$frequency[1:cut];
	cumulative <- cumsum(dsc$frequency)[1:cut];
	
	ephy$marg.probs <- dsc$marg.probs;
 	
 	ephy$shiftnodes <- shiftnodes;
 	ephy$indices <- indices;
 	ephy$frequency <- frequency;
 	ephy$cumulative <- cumulative;
 	ephy$coreshifts <- dsc$coreshifts;
 	ephy$threshold <- threshold;
 	ephy$set.limit <- set.limit;
 	ephy$number.distinct <- length(indices);
 	
	class(ephy) <- 'credibleshiftset';
	return(ephy);	
	
}






