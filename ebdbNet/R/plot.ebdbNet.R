`plot.ebdbNet` <-
function (x, sig.level, interactive = "FALSE", clarify = "TRUE",
	layout = layout.fruchterman.reingold, ...) 
{
	ebdbn <- x
    	if (class(ebdbn) != "ebdbNet") {
        	stop("Error: ", paste(sQuote("ebdbn"), sep = ""), " must be of class ", 
            paste(dQuote("ebdbNet"), sep = ""), sep = "")
    	}

	## Choose edges based on user-defined significance level
	ebdbn.z <- t(ebdbn$z) ## TRANSPOSE FOR VISUALIZATION
	ebdbn.net <- matrix(0, nrow = nrow(ebdbn.z), ncol = ncol(ebdbn.z))
	colnames(ebdbn.net) <- colnames(ebdbn.z)
	rownames(ebdbn.net) <- rownames(ebdbn.z)
	cutoff <- qnorm((1+sig.level)/2)
	ebdbn.net[which(abs(ebdbn.z) > cutoff, arr.ind = TRUE)] <- 1 

	## Feedback graphs
	if(ebdbn$type == "feedback") {
		if(is.null(rownames(ebdbn.net)) == TRUE) rownames(ebdbn.net) <- 1:dim(ebdbn.net)[1];
		if(is.null(rownames(ebdbn.z)) == TRUE) rownames(ebdbn.z) <- 1:dim(ebdbn.z)[1];

		remove.index <- which(colSums(ebdbn.net) == 0 &
			rowSums(ebdbn.net) == 0)
		if(length(remove.index) == dim(ebdbn.net)[1])
			stop("At the current sig.level, no gene-gene interactions have been identified.");
		if(clarify == TRUE & length(remove.index) > 0) 
			ebdbn.net <- ebdbn.net[-remove.index, -remove.index];
		ebdbn.igraph <- graph.adjacency(ebdbn.net, mode=c("directed"), add.rownames = TRUE)
		if(interactive == FALSE) {
			plot(ebdbn.igraph, layout = layout, vertex.label = c(rownames(ebdbn.net)), ...)
		}
		if(interactive == TRUE) {
			tkplot(ebdbn.igraph, vertex.label = c(rownames(ebdbn.net)), ...)
		}
	}

	## Input graphs
	if(ebdbn$type == "input") {
		if(is.null(rownames(ebdbn.net)) == TRUE) rownames(ebdbn.net) <-  paste("u", 1:dim(ebdbn.net)[1], sep = "");
		if(is.null(colnames(ebdbn.net)) == TRUE) colnames(ebdbn.net) <- 1:dim(ebdbn.net)[2];
		if(is.null(rownames(ebdbn.z)) == TRUE) rownames(ebdbn.z) <-  paste("u", 1:dim(ebdbn.net)[1], sep = "");
		if(is.null(colnames(ebdbn.z)) == TRUE) colnames(ebdbn.z) <- 1:dim(ebdbn.net)[2];

		vertex.label.c <- colnames(ebdbn.net)
		vertex.label.r <- rownames(ebdbn.net)
		remove.index.c <- which(colSums(ebdbn.net) == 0)
		remove.index.r <- which(rowSums(ebdbn.net) == 0)
		if(length(remove.index.c) == dim(ebdbn.net)[2])
			stop("At the current sig.level, no input-gene interactions have been identified.");
		if(length(remove.index.r) == dim(ebdbn.net)[1])
			stop("At the current sig.level, no input-gene interactions have been identified.");
		if(clarify == TRUE & length(c(remove.index.c, remove.index.r)) > 0) {
 			if(length(remove.index.r) > 0) {
				ebdbn.net <- ebdbn.net[-remove.index.r,];
				if(is.vector(ebdbn.net) == TRUE)
					ebdbn.net <- matrix(ebdbn.net, nrow = 1);
			}
 			if(length(remove.index.c) > 0) {
				ebdbn.net <- ebdbn.net[,-remove.index.c];
				if(is.vector(ebdbn.net) == TRUE)
					ebdbn.net <- matrix(ebdbn.net, ncol = 1);
			}
			if(length(remove.index.c) > 0) 
				colnames(ebdbn.net) <- colnames(ebdbn.z)[-remove.index.c];
			if(length(remove.index.c) == 0) 
				colnames(ebdbn.net) <- colnames(ebdbn.z);
			if(length(remove.index.r) > 0) 
				rownames(ebdbn.net) <- rownames(ebdbn.z)[-remove.index.r];
			if(length(remove.index.r) == 0) 
				rownames(ebdbn.net) <- rownames(ebdbn.z);
			vertex.label.c <- colnames(ebdbn.net)
			vertex.label.r <- rownames(ebdbn.net)
		}

		tmp <- matrix(which(ebdbn.net != 0, arr.ind = TRUE), ncol = 2)
		tmp[,1] <- tmp[,1] + dim(ebdbn.net)[2]
		if(is.matrix(tmp) == FALSE) tmp <- matrix(tmp, nrow = 1);
		if(is.matrix(tmp) == TRUE) tmp <- as.vector(t(tmp));
		nodes <- c(rep(0, ncol(ebdbn.net)), rep(1, nrow(ebdbn.net)))
		edges <- tmp
		ebdbn.igraph <- graph.bipartite(types = nodes,
			edges = edges, directed = TRUE)

		if(interactive == FALSE) {
			plot(ebdbn.igraph, layout = layout, 
				vertex.color = c(rep("SkyBlue2", ncol(ebdbn.net)),
				rep("green", nrow(ebdbn.net))), 
				vertex.label = c(vertex.label.c, vertex.label.r), ...)
		}
		if(interactive == TRUE) {
			tkplot(ebdbn.igraph, vertex.label = c(vertex.label.c, vertex.label.r), ...)
		}

	}

}

