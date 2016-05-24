prune.cna <- function(x, edges.min=1, size.min=1) {
  
  ##-- Prune nodes based on number of edges and number of members
  ##     prune.cna(net)
  ##
  
  ## Check for presence of igraph package
  oops <- requireNamespace("igraph", quietly = TRUE)
  if (!oops) {
     stop("igraph package missing: Please install, see: ?install.packages")
  }
  
  if(class(x)=="cna") {
    y <- summary.cna(x)
    network=x$community.network
  } else {
    warning("Input should be a 'cna' class object as obtained from cna()")
    network=x
    y <- NULL
  }

  if((edges.min==0) & (size.min==0)){
    stop("Must specify a number greater than 0 for edges.min and/or size.min")
  }
  
  ## Identify nodes with less than 'edges.min' to other nodes.
  nodes.inds <- which(igraph::degree(network) < edges.min)
  
  
  ## Identify nodes with size less than 'size.min'
  ##  cant use V(net$network)$size as these can be scaled
  ##  so we will use the summary information in 'y'
  nodes.inds <- c(nodes.inds, which(y$size < size.min))

  nodes.inds <- unique(nodes.inds)
  
  if( length(nodes.inds) == 0 ) {
    cat( "No Nodes Will Removed based on edges.min and size.min values" )
    output = x
  } else {
    rm.vs <- igraph::V(network)[nodes.inds]
    cat( paste("Removing Nodes:", paste(rm.vs, collapse=", ")),"\n")
    ## Print details of removed with edges
    if(!is.null(y)) {
      w <- cbind(y$tbl[rm.vs,c("id","size")],
                 "edges"=igraph::degree(network)[rm.vs],
                 "members"=y$tbl[rm.vs,c("members")])
      w <- w[order(w$id),]
      write.table(w, row.names=FALSE, col.names=TRUE, quote=FALSE,sep="\t")

      ## Residue raw network
      res2rm <- as.numeric(unlist(y$members[rm.vs]))
      x$communities$membership[res2rm] = NA
    }
  
    d <- igraph::delete.vertices(network, rm.vs)
    
    
    ## Will probably want to keep an edited old community object !!!
    
    output <- list("community.network"=d,
                   "network"= x$network,  ## UNCHANGED!!!
                   "communities"=x$communities)
  }

  class(output) = class(x)
  return(output)
}

