plotNet <-
function(net, fn="", th=1e-6, mylayout=NULL){
  ###require(igraph)
  
  if(is.null(colnames(net))){
    colnames(net) <- paste("V", 1:ncol(net), sep="")
  }
  
  net[abs(net) > th] <- 1
  net[net!=1] <- 0
  diag(net) <- 0
  
  optG <- igraph::graph.adjacency(net,"undirected")
  
  if(is.null(colnames(net))){
    colnames(net) <- paste("V", 1:ncol(net), sep="")
  }
  
  igraph::V(optG)$label <- colnames(net)
  igraph::V(optG)$label.cex <- 10/nrow(net)
  igraph::V(optG)$label.color <- "#060606"
  
  igraph::V(optG)$size <- 500/nrow(net)
  igraph::V(optG)$label.cex = igraph::V(optG)$size*0.04
  igraph::V(optG)$label.font = 2
  
  igraph::V(optG)$frame.color <- NA
  igraph::V(optG)$shape <- "circle"
  igraph::V(optG)$color <- "#0099FF"
  
  igraph::E(optG)$width <- 50/nrow(net) * 2
  igraph::E(optG)$arrow.size <- 0
  igraph::E(optG)$curved <- 0.08
  igraph::E(optG)$color <- "#696A6A"
  
  if(is.null(mylayout)){
    mylayout = igraph::layout.fruchterman.reingold(optG)
  }
  igraph::plot.igraph(optG, layout=mylayout)
  
  if(fn != ""){
	pdf(fn, useDingbats = FALSE)
	igraph::plot.igraph(optG, layout=mylayout)
	dev.off()
	cat(paste("Plot to output file: ", fn, "\n",sep=""))
  }
  
  #title(paste("i=", i, ", Edges=", length(igraph::E(optG)), sep=""))
  return(mylayout)  
}
