ancestors <-
function(y, G, to) { 
  neighbors <- unique(unlist(neighborhood(G, order = vcount(G), nodes = y, mode = "in")))
  v <- V(G)[neighbors]$name
  v <- to[which(to %in% v)]
  return(v)
}
