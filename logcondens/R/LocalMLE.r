LocalMLE <- function(x, w, IsKnot, phi_o, prec){
    n <- length(x)
    res1 <- LocalCoarsen(x, w, IsKnot)
    x2 <- res1$x2
    w2 <- res1$w2
    K <- (1:n) * IsKnot
    K <- K[K > 0]
    res2 <- MLE(x2, w2, phi_o[K])
    phi <- res2$phi
    L <- res2$L
    phi <- LocalExtend(x, IsKnot, x2, phi)
    conv <- as.vector(LocalConvexity(x, phi)) * IsKnot
    Fhat <- LocalF(x, phi)
    H <- 1:n * 0
    JJ <- (1:n) * IsKnot
    JJ <- JJ[JJ > 0]
    
    for (i in 1:(length(JJ) - 1)){
        if (JJ[i + 1] > JJ[i] + 1){
            dtmp <- x[JJ[i + 1]] - x[JJ[i]]
            ind <- (JJ[i] + 1):(JJ[i + 1] - 1)
            mtmp <- length(ind)
            xtmp <- (x[ind] - x[JJ[i]])/dtmp
            wtmp <- w[ind]
            cstmp <- cumsum(xtmp)
            H[ind] <- dtmp * (cumsum(wtmp * xtmp) - xtmp * cumsum(wtmp) + xtmp * sum(wtmp * (1 - xtmp)))
            jtmp1 <- xtmp * J10(phi[ind], phi[JJ[i]] * rep(1, mtmp))
            jtmp2 <- (1 - xtmp) * J10(phi[ind], phi[JJ[i + 1]] * rep(1, mtmp))
            H[ind] <- H[ind] - dtmp^2 * (xtmp * (1 - xtmp)) * (jtmp1 + jtmp2)
        } ## end if
    } ## end for
    
    res <- list(phi = matrix(phi, ncol = 1), L = L, conv = matrix(conv, ncol = 1), H = matrix(H, ncol = 1))
    return(res)
}
