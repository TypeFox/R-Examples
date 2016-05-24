getProb <- function(beta, wn, focusEx){
    focusProb <- exp(beta*wn + focusEx)
    s <- rowSums(focusProb)
    focusProb[which(s == 0 | !is.finite(s)),] <- 1
    focusProb <- focusProb/rowSums(focusProb)
    rMultinom(focusProb)
}

rPotts1Gibbs <- function(nneigh, neighbors, blocks, weights, spatialMat, beta,
                         external, ncolor, colors, newColors){

    for(j in 1:length(blocks)){
        focus <- blocks[[j]]
        focusNei <- neighbors[focus,]
        focusWei <- weights[focus,]
        focusEx <- external[focus,]
        wn <- matrix(0, nrow=nrow(focusNei), ncol=ncolor)
        for(m in 1:nneigh){
            wn <- wn + colors[focusNei[,m],] * focusWei[,m]
        }
        if(!is.null(spatialMat))
            wn <- wn %*% spatialMat
        focusProb <- exp(beta*wn + focusEx)
        s <- rowSums(focusProb)
        focusProb[which(s == 0 | !is.finite(s)),] <- 1
        focusProb <- focusProb/rowSums(focusProb)
        focusCol <- getProb(beta, wn, focusEx)
        newColors[focus] <- focusCol
        colors[focus,] <- do.call(cbind, lapply(1:ncolor, function(x) focusCol==x))
    }
    newColors
}

rPotts1ParDecoup <- function(nvertex, weights, edges, beta,
                             external, ncolor, colors){
    nedge <- nrow(edges)
    bondsTest <- weights *
        ifelse(colors[edges[,1]] - colors[edges[,2]] == 0, 1, 0)
    bondsProb <- 1 - exp(-beta * bondsTest)
    bonds <- edges[runif(nedge) < bondsProb,]
    patches <- getPatches(bonds, nvertex)
    npatch <- length(patches)
    for (j in 1: npatch){
        focus <- patches[[j]]
        if (length(focus)==1)
            focusEx <- external[focus,]
        else
            focusEx <- colSums(external[focus,])
        wn <- t(sapply(1:ncolor, function(k){
            colors[focus] <- k
            sum((1 - weights) * colors[edges[,1]] == colors[edges[,2]])}))
        focusCol <- getProb(beta, wn, focusEx)
        colors[focus] <- focusCol
      
    }

    colors
}


rPotts1 <- function(nvertex, ncolor, neighbors, blocks, edges=NULL,
                    weights=1, spatialMat=NULL, beta, external, colors,
                    algorithm=c("Gibbs", "PartialDecoupling")){
    
    if(nvertex != nrow(external))
        stop("The number of vertex does not match the dimension of 'external'.")

    algorithm <- match.arg(algorithm)
    algorithm <- switch(algorithm, Gibbs = 1, PartialDecoupling = 2)
    
    if(algorithm == 1){
        if(!is.matrix(neighbors))
            stop("'neighbors', a matrix, has to be provided.")
        if(nvertex != nrow(neighbors))
            stop("The number of vertex does not match the dimension of 'neighbors'.")
        dn <- dim(neighbors)
        dw <- dim(weights)
        if(any(dn < dw))
            stop("You provide more 'weights' than 'neighbors'.")
        if(any(dn < dw) || is.null(dw))
            weights <- matrix(as.vector(weights), nrow=dn[1], ncol=dn[2])
    }
    else{
        if(!is.matrix(edges) || ncol(edges) != 2)
            stop("'edges', a matrix of column two, has to be provided.")
        if(nvertex < max(edges))
            stop("The number of vertex has to be at least the same as themax number of vertices in 'edges'.")
        lw <- length(weights)
        le <- nrow(edges)
        if(lw > le)
            stop("You provide more 'weights' than the number of 'edges'.")
        if(lw < le)
            weights <- rep(weights, length=le)
    }

    if(algorithm == 1){
        nneigh <- ncol(neighbors)
        newColors <- rep(0, length(colors))
        colors <- do.call(cbind, lapply(1:ncolor, function(x) colors==x))
        colors <- rbind(colors, rep(0,ncolor))
        rPotts1Gibbs(nneigh, neighbors, blocks, weights, spatialMat, beta,
                     external, ncolor, colors, newColors)
    }
    else{
        rPotts1ParDecoup(nvertex, weights, edges, beta,
                             external, ncolor, colors)
    }

}


