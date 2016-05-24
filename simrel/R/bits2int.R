.bits2int <-
function(x, l2levels){
    nfac <- length(l2levels)
    cumlevels<-c(0,cumsum(l2levels))+1
    intvec <- rep(0,nfac)
    for(i in 1:nfac){
        z <- x[cumlevels[i]:(cumlevels[i+1]-1)]
        intvec[i] <- polyn.eval(z, 2)
    }
    return(intvec)
}
