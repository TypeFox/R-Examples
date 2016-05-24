####################################################################################################
##computes geodesic chains and loops from a predecessor matrix; NOTE: this is an INTERNAL FUNCTION##
##not exported by the package, and as such, it does not provide checks of its argument##############
####################################################################################################
geodesic.chains.loops <-
function(M){
	# M: a predecessor matrix (in Fechnerian scaling context, matrices of the predecessors of the column stimuli
	#    obtained in shortest paths from the row stimuli, as source vertices, to the column stimuli, as target vertices;
	#    using matrices of the psychometric increments as adjacency matrices, of the first and second kind)

	n <- dim(M)[1]
	store <- character(n)
	chains <- matrix(nrow = n, ncol = n)  # geodesic chains
	dimnames(chains) <- dimnames(M)
	diag(chains) <- colnames(M)

	for(i in 1:n){
		store[i] <- chains[i, i]
		done <- NA
		used <- numeric()
		while(any(is.na(chains[i, ]))){
			done <- (which(!is.na(chains[i, ]))[!is.element(which(!is.na(chains[i, ])), used)])[1]
			store[which(is.element(M[i, ], done))] <- paste(store[done], colnames(M)[which(is.element(M[i, ], done))], sep="")
			chains[i, which(is.element(M[i, ], done))] <- store[done]
			used <- append(used, done[!is.element(done, used)])
		}
	}

	loops <- matrix(paste(chains, t(chains), sep=""), ncol = n)  # geodesic loops
	dimnames(loops) <- dimnames(M)
	for(i in 1:n){
		chains[, i] <- paste(chains[, i], colnames(M)[i], sep="")
		loops[i, ] <- paste(loops[i, ], rownames(M)[i], sep="")
	}
	diag(chains) <- diag(loops) <- colnames(M)

	return(list(geodesic.chains = as.data.frame(chains, stringsAsFactors = FALSE),
	            geodesic.loops = as.data.frame(loops, stringsAsFactors = FALSE)))
}
