plot.GMS <-
function(x, fn="", th=1e-6, i=NULL, mylayout=NULL, vars=NULL, ...){
	###require(igraph)
	
	if(is.null(i)){
		cat("Plot optimal network\n")
		optNet <- x$network[[x$opt.index]]
	} else {
		cat("Plot ", i, "-network\n")
		optNet <- x$network[[i]]
	}
	
	optNet[abs(optNet) > th] <- 1
	optNet[optNet!=1] <- 0
	diag(optNet) <- 0
	
	optG <- igraph::graph.adjacency(optNet,"undirected")
	#lout <- layout.fruchterman.reingold(optG)
	#lout <- layout.reingold.tilford(optG)
	
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

	igraph::V(optG)$label <- colnames(optNet)
	igraph::V(optG)$label.cex <- 10/nrow(optNet)
	igraph::V(optG)$label.color <- "#060606"

	igraph::V(optG)$size <- 500/nrow(optNet)
	igraph::V(optG)$label.cex = igraph::V(optG)$size*0.04
	igraph::V(optG)$label.font = 2

	igraph::V(optG)$frame.color <- NA
	igraph::V(optG)$shape <- "circle"
	igraph::V(optG)$color <- "#0099FF"

	igraph::E(optG)$width <- 50/nrow(optNet) * 2
	igraph::E(optG)$arrow.size <- 0
	igraph::E(optG)$curved <- 0.08
	igraph::E(optG)$color <- "#696A6A"

	if(is.null(mylayout)){
		mylayout = igraph::layout.fruchterman.reingold(optG)
	}
	
	if(fn != ""){
		pdf(fn, useDingbats = FALSE)
		#plot(optG, layout=layout.fruchterman.reingold(optG, niter=3000))
		#plot(optG, layout=layout.kamada.kawai(optG, niter=1000))
		igraph::plot.igraph(optG, layout=mylayout)
		#plot(optG,layout=layout.drl)
		dev.off()
		cat(paste("Output file: ", fn, "\n",sep=""))
	}
	if(fn ==""){
		#plot(optG, layout=layout.fruchterman.reingold(optG, niter=3000))
		#plot(optG, layout=layout.kamada.kawai(optG, niter=1000))
		igraph::plot.igraph(optG, layout=mylayout)
	}
	return(mylayout)
}
