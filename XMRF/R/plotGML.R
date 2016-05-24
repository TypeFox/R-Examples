plotGML <-
function(x, fn="", th=1e-6, i=NULL, weight=FALSE, vars=NULL){
	###require(igraph)
	
	if(is.null(i)){
		cat("Write optimal network to GML\n")
		optNet <- x$network[[x$opt.index]]
		i = x$opt.index
	} else {
		cat("Write ", i, "-network to GML\n")
		optNet <- x$network[[i]]
	}
	
	optNet[abs(optNet) > th] <- 1
	optNet[optNet!=1] <- 0
	diag(optNet) <- 0
	
	if(is.null(colnames(optNet))){
		if(is.null(vars)){
			colnames(optNet) <- paste("V", 1:ncol(optNet), sep="")
		}
		else{
			if(length(vars) == ncol(optNet)){
				colnames(optNet) <- vars
			}
			else{
				colnames(optNet) <- paste("V", 1:ncol(optNet), sep="")
			}
		}
	}
	
	G <- igraph::graph.adjacency(optNet,"undirected")
	
	if (weight){
		wmat <- x$D[[i]]
		colnames(wmat) <- colnames(optNet)
		wnet <- wmat * optNet
		G <- igraph::graph.adjacency(wnet, mode="undirected", weighted=TRUE)
	}
	
	igraph::write.graph(G, file=fn , format="gml")
}
