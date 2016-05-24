#############################################################
#
#	branchMeanRateExponential(....)
#
#	Computes time-averaged rate of eponential process

branchMeanRateExponential <- function(t1, t2, p1, p2){
	tol <- 0.00001;
	res <- vector(mode = 'numeric', length = length(t1));
	
	res[which(abs(p2) < tol)] <- p1[which(abs(p2) < tol)];
	nonzero <- which(p2 < -tol);
	p1s <- p1[nonzero];
	p2s <- p2[nonzero];
	t1s <- t1[nonzero];
	t2s <- t2[nonzero];
	res[nonzero] <- (p1s/p2s)*(exp(p2s*t2s) - exp(p2s*t1s)) / (t2s - t1s);
	nonzero <- which(p2 > tol);
	p1s <- p1[nonzero];
	p2s <- p2[nonzero];
	t1s <- t1[nonzero];
	t2s <- t2[nonzero];
	res[nonzero] <- (p1s/p2s)*(2*p2s*(t2s-t1s) + exp(-p2s*t2s) - exp(-p2s*t1s)) / (t2s - t1s);
	return(res);
}