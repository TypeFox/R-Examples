getWeights <- function(mask, neiStruc, format=1){
    if(format==1){
        edges <- getEdges(mask, neiStruc)
        if (is.list(edges))
          edges <- do.call(rbind, edges)
        pos <- matrix(which(mask==1, arr.ind=TRUE), nrow=sum(mask==1))
        weights <- apply(edges, 1, function(x) sqrt(sum((pos[x[1],] - pos[x[2],])^2)))
    }
    else{
        if(format==2){
            neighbors <- getNeighbors(mask, neiStruc)
            nnei <- ncol(neighbors)
            edges <- cbind(rep(1:nrow(neighbors), nnei),
                           as.vector(neighbors))
            pos <- matrix(which(mask==1, arr.ind=TRUE), nrow=sum(mask==1))
            pos <- rbind(pos, rep(NA,ncol(pos)))
            weights <- apply(edges, 1, function(x) sqrt(sum((pos[x[1],] - pos[x[2],])^2)))
            weights[is.na(weights)] <- Inf
            weights <- matrix(weights, ncol=nnei)
        }
        else
            stop("The format has to be 1 or 2")
    }
    1/weights
}
