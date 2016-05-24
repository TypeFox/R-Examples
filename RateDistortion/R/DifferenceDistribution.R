DifferenceDistribution <-
function(channel) {
    ## Compute the difference distribution, P(y - x). This is implemented
    ## only for univariate distributions.
    
    if((ncol(channel$x) > 1) | (ncol(channel$y) > 1)){
        stop("DifferenceDistribution is only defined for univariate distributions.")
    }
    nj <- nrow(channel$x)
    nk <- nrow(channel$y)

    x <- as.vector(channel$x)
    y <- as.vector(channel$y)
    
    z <- sort(unique(as.vector(outer(x, y, FUN = "-"))))
    nz <- length(z)
    pz <- vector(mode = "numeric", length = nz)
    
    eps = 1e-10
    for(k in 1:nz) {

        ## Find all xi and yj such that zk = yj - xi
        for(j in 1:nj) {
            yj <- z[k] + x[j]
            idx <- which(abs(y - yj) < eps)[1]
            if(!is.na(idx)) {
                q <- ConditionalDistribution(channel, j)
                pz[k] <- pz[k] + channel$px[j] * q$p[idx]
            }
        }
    }
    
    list(diff = z, p = pz)
}
