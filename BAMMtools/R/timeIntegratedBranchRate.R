#############################################################
#
#	timeIntegratedBranchRate(....)
#		computes the integral of rates on a branch segment with respect to time
#		Not the average.
#		
#		

timeIntegratedBranchRate <- function(t1, t2, p1, p2){
	tol <- 0.00001;
	res <- vector(mode = 'numeric', length = length(t1));
	# constant rate
	zero <- which(abs(p2) < tol);
	p1s <- p1[zero];
	t1s <- t1[zero];
	t2s <- t2[zero];
	res[zero] <- p1s * (t2s - t1s);
	# declining rate
	nonzero <- which(p2 < -tol);
	p1s <- p1[nonzero];
	p2s <- p2[nonzero];
	t1s <- t1[nonzero];
	t2s <- t2[nonzero];
	res[nonzero] <- (p1s/p2s)*(exp(p2s*t2s) - exp(p2s*t1s));
	# increasing rate
	nonzero <- which(p2 > tol);
	p1s <- p1[nonzero];
	p2s <- p2[nonzero];
	t1s <- t1[nonzero];
	t2s <- t2[nonzero];
	res[nonzero] <- (p1s/p2s)*(2*p2s*(t2s-t1s) + exp(-p2s*t2s) - exp(-p2s*t1s));
	return(res);
}
