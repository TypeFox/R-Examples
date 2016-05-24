laplacian <-
function (size = 16, compose = FALSE) 
{
    gmat <- matrix(0, size, size)
    xx <- seq(size)
    for (v in xx) gmat[, v] <- sqrt(2/size) * cos(((xx - 0.5) * 
        pi * (v - 1))/size)
    gmat[, 1] <- gmat[, 1]/sqrt(2)
    lvec <- -(2 * size^2) * (1 - cos(((xx - 1) * pi)/size))
    gmat <- kronecker(gmat, gmat)
    lvec <- rep(lvec, rep(size, size)) + rep(lvec, size)
    if (compose) 
        gmat %*% (lvec^2 * t(gmat))
    else list(vectors = gmat, values = lvec^2)
}

