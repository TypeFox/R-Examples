sc.components <-
function(D, to) {
  from <- NULL
  A <- as.matrix(get.adjacency(D))
  v <- get.vertex.attribute(D, "name")
  s <- v[which(vertex.attributes(D)$description == "S")]
  e <- E(D)
  bidirected <- c()
  selection <- E(D)[from(s)]
  indices <- which(A >= 1 & t(A) >= 1, arr.ind = TRUE)
  sc.comp <- list()
  if (nrow(indices) > 0) {
    for (i in 1:nrow(indices)) {
      bidirected <- c(bidirected, e[v[indices[i,1]] %->% v[indices[i,2]]])
    }
  }
  D.bidirected <- subgraph.edges(D, union(bidirected, selection), delete.vertices = FALSE)
  subgraphs <- decompose.graph(D.bidirected)
  for (i in 1:length(subgraphs)) { 
    v.sub <- get.vertex.attribute(subgraphs[[i]], "name")
    v.sub <- to[which(to %in% v.sub)]
    sc.comp[[i]] <- v.sub 
  }
  return(sc.comp)
}
