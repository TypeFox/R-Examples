`ecol.dist` <-
function(x, method=sorenson, type="dis") {
    N<-ncol(x)
    mat.names <- list(colnames(x),colnames(x))
    sim<- matrix(,N,N,dimnames=mat.names)
    for (i in 2:N) {
        for (j in 1:(i-1)) sim[i,j] <- method(x[,i],x[,j])
    }
    if (type=="sim") typ <- "similarity" 
    else {
        typ <- "dissimilarity"
        sim<-1-sim
    }
    sim.mat<-as.dist(sim)
    attr(sim.mat, "Size") <- N
    attr(sim.mat, "call") <- match.call()
    attr(sim.mat, "type") <- typ
    return(sim.mat)
}

