#############################################################
#
#	getCohortMatrix(....)
#
# 	Each entry of this matrix represents the expected 
#	probability that a pair[i, j] of tips will have the same 
#	rate parameters due to BAMM model.
#
# 	Should modify this to allow exponential, spherical,
#		and other possible correlation structures.
#	Need to make a corStruct class that works with this
#		for GLS analyses

# getCohortMatrix <- function(ephy) {

	# if (!'bammdata' %in% class(ephy)) {
		# stop("Object ephy must be of class bammdata\n");
	# }
	
	# TOL <- 0.0001;
	# corMat <- matrix(0, nrow=length(ephy$tip.label), ncol=length(ephy$tip.label));
	# n <- length(ephy$numberEvents);
	# for (i in 1:length(ephy$tipStates)) {
		# dd <- dist(ephy$tipStates[[i]]);
		# cmat <- as.matrix(dd);	
		# corMat <- corMat + (cmat < TOL)/n;
		
	# }
	# #rownames(corMat) <- ephy$tip.label;
	# #colnames(corMat) <- ephy$tip.label;	
	# dimnames(corMat)[1:2] <- list(ephy$tip.label);
	# #return(corMat/length(ephy$numberEvents));
	# return(corMat);
# }

getCohortMatrix <- function(ephy) {
	if (!'bammdata' %in% class(ephy)) {
		stop("Object ephy must be of class bammdata\n");
	}
	tipStates <- unlist(ephy$tipStates);
	Ntips <- length(ephy$tip.label);
	Nsamples <- length(ephy$tipStates);
	mat <- .C("cohort_matrix",as.integer(tipStates), as.integer(Nsamples), as.integer(Ntips), double(Ntips*Ntips), PACKAGE = "BAMMtools")[[4]];
	dim(mat) <- c(Ntips, Ntips);
	dimnames(mat) <- list(ephy$tip.label, ephy$tip.label);
	diag(mat) <- rep(1.0,Ntips);
	return(mat);
}