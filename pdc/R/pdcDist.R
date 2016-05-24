pdc.dist <- function(X, m=NULL, t=NULL, divergence=symmetricAlphaDivergence)
{
	warning("Call to pdc.dist(...) is deprecated!")
	return(pdcDist(X, m, t, divergence));
}

pdcDist <-
function(X, m=NULL, t=NULL, divergence=symmetricAlphaDivergence)
{
	if (is.null(t) | is.null(m)) {
		ent <- entropyHeuristic(X)
		
		if (is.null(m)) {
			m <- ent$m;
		}
		
		if (is.null(t)) {
			t <- ent$t;
		}
	}
	

	if (length(dim(X)) == 2) { 
		
		codebooks <- convertMatrix(X,m,t);
		D <- divergenceMatrix( codebooks, divergence );
	
	} else if (length(dim(X))==3) {

		codebooks <- convertMatrixMultichannel(X,m,t);
		num.channels <- dim(X)[3]
		D <- divergenceMatrixMultichannel( codebooks, divergence, num.channels );		
		
	} else {
		stop("Invalid dimensionality of data object!");
	}
	
	pdcdist <- as.dist(D);
	
	return(pdcdist);
}
