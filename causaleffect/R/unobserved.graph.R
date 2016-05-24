unobserved.graph <-
function(G) {
  unobs.edges <- which(edge.attributes(G)$description == "U")
  u <- length(unobs.edges)
  if (u > 0) {
    new.nodes <- paste0("U", 1:u)
    G <- G + vertices(new.nodes)
    v <- get.vertex.attribute(G, "name")
    e <- get.edges(G, unobs.edges)
    for (i in 1:u) G <- G + edges(c(new.nodes[i], v[e[i,1]]), c(new.nodes[i], v[e[i,2]]))
    unobs.edges <- which(edge.attributes(G)$description == "U")
    obs.edges <- setdiff(E(G), E(G)[which(edge.attributes(G)$description == "U")])
    G.unobs <- subgraph.edges(G, E(G)[obs.edges], delete.vertices = FALSE)
    return(G.unobs)
  }
  return(G)
}
