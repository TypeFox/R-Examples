BlocksGibbs <- function(n, nvertex, ncolor, neighbors, blocks,
                        weights=1, spatialMat=NULL, beta){
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

    nneigh <- ncol(neighbors)
    oneIteration <- sample(1:ncolor, size=nvertex, replace=TRUE)
    oneIteration <- do.call(cbind, lapply(1:ncolor, function(x) oneIteration==x))
    oneIteration <- rbind(oneIteration, rep(0,ncolor))
    colors <- matrix(0, nrow=nvertex, ncol=n)
    for(i in 1:n)
        for(j in 1:length(blocks)){
            focus <- blocks[[j]]
            focusNei <- neighbors[focus,]
            focusWei <- weights[focus,]
            wn <- matrix(0, nrow=nrow(focusNei), ncol=ncolor)
            for(m in 1:nneigh){
                wn <- wn + oneIteration[focusNei[,m],] * focusWei[,m]
            }
            if(!is.null(spatialMat))
                wn <- wn %*% spatialMat
            focusProb <- exp(beta * wn)
            s <- rowSums(focusProb)
            focusProb[which(s == 0 | !is.finite(s)),] <- 1
            focusProb <- focusProb / rowSums(focusProb)
            focusCol <- rMultinom(focusProb)
            colors[,i][focus] <- focusCol
            oneIteration[focus,] <- do.call(cbind, lapply(1:ncolor, function(x) focusCol==x))
        }
    colors
}


