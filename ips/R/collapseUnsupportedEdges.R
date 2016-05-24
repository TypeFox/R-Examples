collapseUnsupportedEdges <- function(phy, value, cutoff){
  
  if ( missing(value) ) value <- "node.label"
  
  stat <- as.numeric(phy[[value]])
  collapse <- which(stat < cutoff) + Ntip(phy)
  ## collapse nodes in post-order traversal!!
  ## ----------------------------------------
  for ( i in rev(collapse) ){
    id <- phy$edge[, 2] == i
    id2 <- phy$edge[, 1] == i
    phy$edge[id2, 1] <- phy$edge[id, 1]
    phy$edge <- phy$edge[!id, ]
    phy$edge[phy$edge > i] <- phy$edge[phy$edge > i] - 1
    
    ## edge lengths
    phy$edge.length[id2] <- phy$edge.length[id2] + phy$edge.length[id]
    phy$edge.length <- phy$edge.length[!id]
    phy$node.label <- phy$node.label[!id]
    phy$Nnode <- phy$Nnode - 1
    #phy <- fixNodes(phy)
  }
  phy
}