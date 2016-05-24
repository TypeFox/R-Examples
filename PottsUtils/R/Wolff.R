findBonds<- function(nei, neiCol, col, neiWei, beta){
    bondsTest  <- ifelse(neiCol - col == 0, 1, 0)
    bondsProb <- 1 - exp(-beta * neiWei * bondsTest)
    nei[runif(length(nei)) < bondsProb]
}

getCluster <- function(oneIteration, neighbors, weights, nvertex, beta){
    pick <- sample(1:nvertex, 1)
    cluster <- pick
    queue <- pick

    while(length(queue) > 0){
        first <- queue[1]
        queue <- queue[-1]
        firstNei <- neighbors[first,]
        firstWei <- weights[first,]
        firstNeiCol<- oneIteration[firstNei]
        firstCol <- oneIteration[first]
        bonds <- findBonds(firstNei, firstNeiCol, firstCol, firstWei, beta)
        if(length(bonds) > 0){
            new <- bonds[is.na(match(bonds, cluster))]
            cluster <- c(cluster, new)
            queue <- c(queue, new)
        }
    }

    cluster
}

Wolff <- function(n, nvertex, ncolor, neighbors, weights=1, beta){
    if(nvertex != nrow(neighbors))
        stop("The number of vertex does not match the dimension of 'neighbors'.")
    if (!is.matrix(neighbors))
        stop("'neighbors' has to be a matrix.")
    dn <- dim(neighbors)
    dw <- dim(weights)
    if (any(dn < dw))
        stop("You provide more 'weights' than 'neighbors'.")
    if (any(dn < dw) || is.null(dw))
        weights <- matrix(as.vector(weights), nrow=dn[1], ncol=dn[2])

    colors <- matrix(0, nrow=nvertex+1, ncol=n)
    oneIteration <- c(sample(x=1:ncolor, nvertex, replace=TRUE), ncolor+1)
    for(i in 1:n){
        cluster <- getCluster(oneIteration, neighbors, weights, nvertex, beta)
        kick <- oneIteration[cluster[1]]
        if (ncolor==2)
            newColor <- (1:ncolor)[-kick]
        else
            newColor <- sample((1:ncolor)[-kick],1)
        oneIteration[cluster] <- newColor
        colors[,i] <- oneIteration
    }

    colors[-(nvertex+1),]

}
