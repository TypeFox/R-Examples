`mNNtest` <-
function (info, obsN = NULL) 
{
    if (is.null(obsN)) 
        obsN <- as.vector(t(info$ON))
    expN <- as.vector(t(info$EN))
    varN <- diag(info$VarN)
    Z <- (obsN - expN)/sqrt(varN)
    names(Z) <- dimnames(info$VarN)[[1]]
    k <- nrow(info$EN)
    delta <- as.matrix(obsN - expN)
    C <- t(delta) %*% ginv(info$VarN) %*% delta
    pC <- 1 - pchisq(C, k * (k - 1))
    Ci <- rep(0, k)
    for (i in 1:k) {
        i1 <- 1 + (i - 1) * k
        i2 <- i1 + (k - 2)
        Ci[i] <- t(delta[i1:i2, ]) %*% solve(info$VarN[i1:i2, 
            i1:i2]) %*% delta[i1:i2]
    }
    pCi <- 1 - pchisq(Ci, k - 1)
    list(Z = Z, C = c(C, pC), Ci = cbind(Ci, pCi))
}

