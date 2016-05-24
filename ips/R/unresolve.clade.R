unresolve.clade <- function(phy, node){
  tip <- descendants(phy, node, type = "terminal")
  del <- vector()
  for ( i in tip ){
    # assign sum of edge lengths to terminal edge
    n <- i
    e <- vector()
    while ( n[1] != node ){
      e <- c(which.edge(phy, n[1]), e)
      n <- c(phy$edge[phy$edge[, 2] == n[1], 1], n)
    }
    if( !is.null(phy$edge.length) )
      phy$edge.length[tail(e, 1)] <- sum(phy$edge.length[e])
    del <- c(del, head(e, -1))
  }
  phy$edge[phy$edge[, 2] %in% tip, 1] <- node
  del <- unique(del)
  # adjust node index
  id <- phy$edge[del, 2]
  id <- which(phy$edge > id, arr.ind = TRUE)
  phy$edge[id] <- phy$edge[id] - length(del)
  
  phy$edge <- phy$edge[-del, ]
  if( !is.null(phy$edge.length) )
    phy$edge.length <- phy$edge.length[-del]
  phy$Nnode <- nrow(phy$edge) + 1 - Ntip(phy) 
  phy
}