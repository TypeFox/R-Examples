`K012` <-
function (X, fijo, i, j, nsim = 99, nrank = 1, r = NULL, correction = "isotropic") 
{
    marx <- marks(X)
    fijo <- (marx == fijo)
    I <- (marx == i)
    J <- (marx == j)
    cosa <- Kmulti.ls(X, fijo, I, r, corre = correction)
    r = cosa$r
    k01.obs <- cosa[[3]]
    k02.obs <- Kmulti.ls(X, fijo, J, r, corre = correction)[[3]]
    k01.sim <- NULL
    k02.sim <- NULL
    cat("Generating simulations...")
    for (n in 1:nsim) {
        progressreport(n, nsim)
        X$marks[I | J] <- sample(X$marks[I | J])
        marx <- marks(X)
        I <- (marx == i)
        J <- (marx == j)
        k01.sim <- cbind(k01.sim, Kmulti.ls(X, fijo, I, r, corre = correction)[[3]])
        k02.sim <- cbind(k02.sim, Kmulti.ls(X, fijo, J, r, corre = correction)[[3]])
    }
    orderstat <- function(x, n) sort(x)[n]
    k01.lo <- apply(k01.sim, 1, orderstat, n = nrank)
    k01.hi <- apply(k01.sim, 1, orderstat, n = nsim - nrank + 
        1)
    k02.lo <- apply(k02.sim, 1, orderstat, n = nrank)
    k02.hi <- apply(k02.sim, 1, orderstat, n = nsim - nrank + 
        1)
    k01 <- cosa
    k01[[2]] <- k01.hi
    k01[[3]] <- k01.obs
    k01[[4]] <- k01.lo
    attributes(k01)$labl <- c(attributes(cosa)$labl[1], "hi(r)", 
        attributes(cosa)$labl[3], "lo(r)")
    attributes(k01)$names <- c(attributes(cosa)$names[1], "hi", 
        attributes(cosa)$names[3], "lo")
    attributes(k01)$names <- c(attributes(cosa)$names[1], "hi", 
        attributes(cosa)$names[3], "lo")
    attributes(k01)$desc <- c(attributes(cosa)$desc[1], "upper pointwise envelope of simulations", 
        attributes(cosa)$desc[3], "lower pointwise envelope of simulations")
    k02 <- cosa
    k02[[2]] <- k02.hi
    k02[[3]] <- k02.obs
    k02[[4]] <- k02.lo
    attributes(k02)$labl <- c(attributes(cosa)$labl[1], "hi(r)", 
        attributes(cosa)$labl[3], "lo(r)")
    attributes(k02)$names <- c(attributes(cosa)$names[1], "hi", 
        attributes(cosa)$names[3], "lo")
    attributes(k02)$names <- c(attributes(cosa)$names[1], "hi", 
        attributes(cosa)$names[3], "lo")
    attributes(k02)$desc <- c(attributes(cosa)$desc[1], "upper pointwise envelope of simulations", 
        attributes(cosa)$desc[3], "lower pointwise envelope of simulations")
    return(list(k01 = k01, k02 = k02))
}

