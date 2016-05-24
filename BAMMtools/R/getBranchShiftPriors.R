
getBranchShiftPriors <- function(phy, expectedNumberOfShifts) {
	
	Nmax <- 1000;
	
	geom_p <- 1 / (expectedNumberOfShifts + 1);
	prior <- dgeom(1:Nmax, geom_p);
 
	pvec <- phy$edge.length / sum(phy$edge.length);
	
	pp <- numeric(length(phy$edge.length));
 
	for (i in 1:length(prior)){
		# probability of getting 0 shifts on branch given ns total 
		#  given the branch lengths etc
		#	weighted by the probability of that shift category
		
		pp <- pp + (1 - dbinom(0, i, prob=pvec)) * prior[i];
		
	}
	
	obj <- phy;	
	obj$edge.length <- pp;
	return(obj);
}

