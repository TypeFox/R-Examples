getEdges <- function(mask, neiStruc){
    neighbors <- getNeighbors(mask, neiStruc)
    nvertex <- nrow(neighbors)
    nneighbor <- ncol(neighbors)
    edges <- do.call(rbind, lapply(1:nneighbor, function(i){
        edges <- cbind(1:nvertex, neighbors[,i])
        edges[edges[,1] < edges[,2],]
        }))
    edges <- edges[edges[,2] != nvertex+1, ]
    edges
}
