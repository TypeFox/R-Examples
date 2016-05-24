mritc.icm <- function(y, neighbors, blocks, spatialMat=diag(1,3), beta=0.4, mu, sigma, err=1e-4, maxit=20, verbose){

    checkErrors(mu=mu, sigma=sigma, err=err)

    if (length(err) < 2) err <- rep(err, length.out = 2)

    k <- length(mu)
    nvert <- length(y)
    yunique <- sort(unique(y))
    n.yunique <- length(yunique)
    ymatch <- match(y, yunique)
    nvert <- length(y)
    nneigh <- ncol(neighbors)

    check <- getCheck(nneigh, k, beta, spatialMat)
    indices <- initialIndices(y, nvert, mu, sigma, k, sub=FALSE)

    prob <- matrix(0, ncol=k, nrow=nvert)
    it <- 0

    repeat{
        it <- it + 1
        #update indices
        muold <- mu
        sigmaold <- sigma
        den <- getDen(yunique, n.yunique, ymatch, mu, sigma)

        
        for (i in 1:length(blocks)){
            points <- blocks[[i]]
            prob[points,] <- updateProbs(neighbors[points,], nneigh, indices,
                                 den[points,], k, check)
            indices.focus <- max.col(prob[points,])
            indices.focus <- do.call(rbind, lapply(1:k, function(x) indices.focus==x))
            indices[points,] <- t(indices.focus)
        }

        #update mu and sigma
        Nj <- colSums(indices)
        mu <- colSums(y*(indices[-(nvert+1),]))/Nj
        sigma <- sqrt(colSums(
                     (y*(indices[-(nvert+1),])- t(t(indices[-(nvert+1),])*mu))^2)/Nj)

        flag <- checkStopVerbose(muold, mu, sigmaold, sigma, err, it, maxit, verbose)
        if(flag==1) break

    }

    list(prob=prob/rowSums(prob), mu=mu, sigma=sigma)
}
