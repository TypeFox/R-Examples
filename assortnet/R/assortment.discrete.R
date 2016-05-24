
assortment.discrete <- function(graph, types, weighted=TRUE, SE=FALSE, M=1) {
	
	graph <- as.matrix(graph)
	
	if (!weighted) {
		graph[graph > 0] <- 1
		graph[graph < 1] <- 0
	}
	
	total_types <- unique(types)
	
	out <- matrix(0,nrow=length(total_types),ncol=length(total_types))

	for (i in 1:length(total_types)) {
		for (j in 1:length(total_types)) {
			out[i,j] <- sum(graph[which(types==total_types[i]),which(types==total_types[j])])/sum(graph)
		}
	}
	
	r <- ( sum(diag(out)) - sum(rowSums(out)*colSums(out)) ) / (1 - sum(rowSums(out)*colSums(out)) )
	

	#if (SE) {
	#	
	#	se <- (sum(rowSums(out)*colSums(out)) + sum(rowSums(out)*colSums(out))^2 - sum((rowSums(out)^2)*colSums(out)) - sum(rowSums(out)*(colSums(out)^2))) / (1 - sum(rowSums(out)*colSums(out)) )
	#	se <- se/sum(graph>0)
	#
	#}
	
	if (SE) {	
		se <- 0
		
		N <- which(graph>0)
		E <- seq(1,length(N),M)
		if (E[length(E)] < length(N))
			E <- c(E,(length(N)+1))
		
		for (g in 1:(length(E)-1)) {
			graph2 <- graph
			graph2[N[E[g]:(E[g+1]-1)]] <- 0
	
			out2 <- matrix(0,nrow=length(total_types),ncol=length(total_types))

			for (i in 1:length(total_types)) {
				for (j in 1:length(total_types)) {
					out2[i,j] <- sum(graph2[which(types==total_types[i]),which(types==total_types[j])])/sum(graph2)
				}
			}
	
			ri <- ( sum(diag(out2)) - sum(rowSums(out2)*colSums(out2)) ) / (1 - sum(rowSums(out2)*colSums(out2)) )
		
			se <- se + (ri-r)^2
		}
		
		se <- sqrt(((length(E)-1)/length(E))*se)
		
	}
	
	out <- cbind(out,rowSums(out))
	out <- rbind(out,colSums(out))
	colnames(out) <- c(as.character(total_types),"ai")
	rownames(out) <- c(as.character(total_types),"bi")
	
	if (SE) {	
		return(list(r=r,se=se,mixing_matrix=out))
	} else {
		return(list(r=r,mixing_matrix=out))
	}
	
}
