c.components <-
function(G, to) {
  A <- as.matrix(get.adjacency(G))
  v <- get.vertex.attribute(G, "name")
  e <- E(G)
  bidirected <- c()
  indices <- which(A >= 1 & t(A) >= 1, arr.ind = TRUE)
  c.comp <- list()
  if (nrow(indices) > 0) {
    for (i in 1:nrow(indices)) {
      bidirected <- c(bidirected, e[v[indices[i,1]] %->% v[indices[i,2]]])
    }
  }
  G.bidirected <- subgraph.edges(G, bidirected, delete.vertices = FALSE)
  subgraphs <- decompose.graph(G.bidirected)
  for (i in 1:length(subgraphs)) { 
    v.sub <- get.vertex.attribute(subgraphs[[i]], "name")
    v.sub <- to[which(to %in% v.sub)]
    c.comp[[i]] <- v.sub 
  }
  return(c.comp)
}
