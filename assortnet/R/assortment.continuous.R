
assortment.continuous <- function(graph, vertex_values, weighted=TRUE, SE=FALSE, M=1) {

	graph <- as.matrix(graph)
	
	if (!weighted) {
		graph[graph > 0] <- 1
		graph[graph < 1] <- 0
	}
	
	i  <- which(graph>0,arr.ind=TRUE)
	ji <- which(graph>0,arr.ind=TRUE)[,1]
	ki <- which(graph>0,arr.ind=TRUE)[,2]
	top1 <- sum(graph[i]*vertex_values[ji]*vertex_values[ki])
	top2 <- sum(graph[i]*vertex_values[ji])
	top3 <- sum(graph[i]*vertex_values[ki])
	bottom1 <- sum(graph[i]*vertex_values[ji]^2)
	bottom2 <- sum(graph[i]*vertex_values[ji])^2
	bottom3 <- sum(graph[i]*vertex_values[ki]^2)
	bottom4 <- sum(graph[i]*vertex_values[ki])^2
	
	r<- (top1 - (1/sum(graph) * top2 * top3)) / sqrt((bottom1-(1/sum(graph))*bottom2)*(bottom3-(1/sum(graph))*bottom4))

	if (SE) {	
		se <- 0
		
		N <- which(graph>0)
		E <- seq(1,length(N),M)
		if (E[length(E)] < length(N))
			E <- c(E,(length(N)+1))
		
		for (g in 1:(length(E)-1)) {
			graph2 <- graph
			graph2[N[E[g]:(E[g+1]-1)]] <- 0
	
			i  <- which(graph2>0,arr.ind=TRUE)
			ji <- which(graph2>0,arr.ind=TRUE)[,1]
			ki <- which(graph2>0,arr.ind=TRUE)[,2]
			top1 <- sum(graph2[i]*vertex_values[ji]*vertex_values[ki])
			top2 <- sum(graph2[i]*vertex_values[ji])
			top3 <- sum(graph2[i]*vertex_values[ki])
			bottom1 <- sum(graph2[i]*vertex_values[ji]^2)
			bottom2 <- sum(graph2[i]*vertex_values[ji])^2
			bottom3 <- sum(graph2[i]*vertex_values[ki]^2)
			bottom4 <- sum(graph2[i]*vertex_values[ki])^2

			ri<- (top1 - (1/sum(graph2) * top2 * top3)) / sqrt((bottom1-(1/sum(graph2))*bottom2)*(bottom3-(1/sum(graph2))*bottom4))
		
			se <- se + (ri-r)^2
		}
		
		se <- sqrt(((length(E)-1)/length(E))*se)
		
		return(list(r=r, se=se))
	} else {
		return(list(r=r))
	}

}
