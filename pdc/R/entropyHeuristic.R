entropy.heuristic <-
function(X, m.min=3, m.max=7, t.min=1, t.max=1)
{
	warning("Call to entropy.heuristic(...) is deprecated!")
	return(entropyHeuristic(X, m.min, m.max, t.min, t.max))
}

entropyHeuristic <-
function(X, m.min=3, m.max=7, t.min=1, t.max=1)
{
	X <- as.matrix(X)

	if (m.max>7) {
	 # issue a message
	 message(paste("No fast implementation for embedding sizes larger than 7 is available. Calculation might be slow!"));	
	}

	ent <- matrix(rep(0, (m.max-m.min+1)*(t.max-t.min+1)*3 ),ncol=3 )
	k<-1
	for (j in t.min:t.max)
	 for (i in m.min:m.max)
	 {
	 	ent[k,1] <- j
	 	ent[k,2] <- i
		ent[k,3] <- mean(apply(FUN=codebook.entropy, MARGIN=2, X, m=i,t=j))
		k <- k+1
	 }
	best <- which.min(ent[,3]);
	
	result <- list()
	result$entropy.values <- ent
	result$m <- ent[best,2]
	result$t <- ent[best,1]
	result$m.range <- m.min:m.max
	result$t.range <- t.min:t.max
	
	class(result) <- "mine"
	
	return (result);
}
