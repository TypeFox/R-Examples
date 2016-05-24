TopoS <- function(EV.dist, Cells.dist){
    dim(Cells.dist) <- NULL
    dim(EV.dist) <- NULL
    beta.est <- (1/(Cells.dist %*% Cells.dist)) %*% Cells.dist %*% EV.dist
    Cells.dist.est <- (Cells.dist * beta.est) - EV.dist
    return(1 - sqrt((Cells.dist.est %*% Cells.dist.est) / (EV.dist %*% EV.dist)))
}

distmirr <- function(dis){
    n <- attr(dis, "Size")
    mirr <- matrix(0, n, n)
    mirr[lower.tri(mirr)] <- dis
    mirr + t(mirr)
}
