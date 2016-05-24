comp2 <-
function (y, z) 
{
    if (length(y) != length(z)) 
        stop("Vectors y and z should have equal length!")
    n <- length(y)
    matY <- outer(y, y, "==")
    matZ <- outer(z, z, "==")
    mat <- matY & matZ
    jstar <- led <- rep(FALSE, n)
    for (j in seq(n)) {
        if (!led[j]) {
            jstar[j] <- TRUE
            if (j == n) 
                break
            ji <- seq(j + 1, n)
            tt <- mat[ji, j] %in% TRUE
            led[ji] <- led[ji] | tt
        }
        if (all(led[-seq(j)])) 
            break
    }
    ystar <- y[jstar]
    zstar <- z[jstar]
    nstar <- apply(as.matrix(mat[, jstar]), 2, sum)
    rstar <- length(nstar)
    idx <- match(y, ystar)
    return(list(ystar = ystar, zstar = zstar, nstar = nstar, 
        rstar = rstar, idx = idx))
}
