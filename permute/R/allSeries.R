`allSeries` <- function(n, nperms, mirror = FALSE)
{
    v <- seq_len(n)
    X <- matrix(nrow = nperms, ncol = n)
    for(i in v) {
        X[i,] <- seq(i, length = n)%%n + 1
    }
    ## if mirroring, rev the cols of X[v,]
    ## but only if nperms > 2
    if(mirror && (nperms > 2))
        X[(n+1):(2*n),] <- X[v, rev(v)]
    X
}
