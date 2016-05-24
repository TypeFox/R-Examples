###############################################################################
## convert similarity matrix to dist object
###############################################################################
sim2dist <- function(x, maxSim = 1){
    stopifnot(is.matrix(x))
    stopifnot(isSymmetric(x))
    stopifnot(maxSim >= max(x))
    d <- maxSim - x # from similarity to distance
    return(as.dist(d))
}
