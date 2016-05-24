#Update probs
updateProbs <- function(neighbors, nneigh, indices, den, k, check){
    Nij <- matrix(0, nrow=nrow(neighbors), ncol=k)
    for(m in 1:nneigh){
        Nij <- Nij + indices[neighbors[,m],]
    }
    or <- Nij[,1]
    for(i in 2:k)
        or <- or + Nij[,i]*(nneigh+1)^(i-1)
    or <- or + 1
    prob <- den*check[or,]
    prob <- prob/rowSums(prob)
    prob[is.na(rowSums(prob)),] <- rep(1/k, k)
    prob[prob==Inf] <- 1/k
    prob
}

#Update mus
updateMus <- function(prob, y){
    inten <- y * prob
    mu <- colSums(inten)/colSums(prob)
    mu
}

#Update sds
updateSds <- function(prob, y, mu, nvert, k){
    mu <- matrix(mu, nrow=nvert, ncol=k, byrow=T)
    diff2 <- (y-mu)^2*prob
    sigma2 <- colSums(diff2)/colSums(prob)
    sqrt(sigma2)
}

updateIndicesHMRFEM <- function(neighbors, nneigh, blocks, nblocks, k, indices, check, den){
    .Call("updateIndicesHMRFEM", blocks, neighbors, nneigh, k, indices, check, den)
}

