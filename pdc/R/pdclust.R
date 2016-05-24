pdclust <-
function(X, m=NULL, t=NULL, divergence=symmetricAlphaDivergence, clustering.method="complete") 
{
	user.m <- !is.null(m);
	user.t <- !is.null(t);
	
	if ((is.null(m)) && (is.null(t))) {
			m <- entropyHeuristic(X)$m;
	}
	
	if (is.null(m)) {
			m <- entropyHeuristic(X, t.min=t, t.max=t)$m;
	}
	
	if (is.null(t)) {
			t <- 1;
	}
	


	# calculate divergence matrix			
	D <- pdcDist(X,m,t,divergence);
	
	# start hierarchical clustering
	if (clustering.method == "complete") {
		hcl <- hclust(D, method="complete")
	} else if (clustering.method == "average") {
		hcl <- hclust(D, method="average")
	}  else if (clustering.method == "single") {
		hcl <- hclust(D, method="single")
	} else {
		stop("Invalid clustering method!")
	}
	
	# add meta info
	hcl$divergence <- divergence
	hcl$m <- m
	hcl$t <- t
	hcl$user.specified.m <- user.m
	hcl$user.specified.t <- user.t
	hcl$N <- length(hcl$order)
	hcl$data <- X
	hcl$D <- D
	hcl$multichannel <- length(dim(X))==3
	hcl$labels <- colnames(X)
	
	# wrap hclust result
	class(hcl) <- c("pdclust","hclust")
	
	return(hcl);	
}
