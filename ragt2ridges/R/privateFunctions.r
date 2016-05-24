.is.int <- function (x, tolerance = .Machine$double.eps){
    abs(x - round(x)) < tolerance
}


.LL <- function (S, P){
    LL <- -log(det(P)) + .trace(S %*% P)
    return(LL)
}

.trace <- function (M) {
    return(sum(diag(M)))
}


.cvlPchordal <- function(lambda, Y,  target=default.target(covML(Y)),  zeros, cliques, separators, type="Alt"){
	slh <- numeric()
	for (i in 1:nrow(Y)){
		S <- covML(Y[-i, ])
		slh <- c(slh, .LL(crossprod(Y[i, , drop=FALSE]), ridgePchordal(S, lambda, target=target, zeros=zeros, cliques=cliques, separators=separators, type=type, verbose=FALSE)))
	}
	return(mean(slh))        
}



.zeros2chordalCIG <- function(zeros, nCovariates){
	####################################################################
	# converts CI pattern of error precision (specified as zeros)
	# into graph-object, which is triangulated and decomposed
	####################################################################

	# zeros to adjacency matrix
	adjMat <- matrix(1, nCovariates, nCovariates)
	adjMat[zeros] <- 0
	adjMat[cbind(zeros[,2], zeros[,1])] <- 0	
	diag(adjMat) <- 0

    	# convert adjacency into graphNel object
	G <- graph.adjacency(adjMat, mode="undirected")
	G <- igraph.to.graphNEL(G)

	# is graph chordal?
	if (!is.triangulated(G)){ G <- triangulate(G) }

        # decompose in cliques and separators
        decomp <- rip(G)
   
	# return decomposition
	return(decomp)
}


.isStationary <- function(A){
	#############################################################################
	# tests whether the VAR(1) process associated with matrix A is stationary
	#############################################################################
	
	# assess stationarity	
    	evs <- abs(eigen(A, only.values=TRUE)$values)
    	if (max(evs) > 1){ print("WARNING: the VAR(1) process associated with the supplied matrix A is not stable!") }
}

