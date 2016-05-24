HessSurvPH <-
function (WW, Y, Y2, ew) {
    WY <- if (is.null(WW)) Y2 else cbind(WW[indT, , drop = FALSE], Y2)
    ncwy <- ncol(WY)
    p <- ncww + 1
    H <- matrix(0, p, p)
    index <- which(lower.tri(H, TRUE), arr.ind = TRUE)
    nn <- nrow(index)
    h <- numeric(nn)
    for (i in 1:nn) {
        i1 <- index[i, 1]
        if (i1 == p)
            i1 <- seq(p, ncwy)
        i2 <- index[i, 2]
        if (i2 == p)
            i2 <- seq(p, ncwy)
        pp <- rowsum(lambda0. * (WY[, i1] * WY[, i2]) * ew, indT)
        h[i] <- - sum((pp * p.byt.) %*% wGH)
    }
    H[lower.tri(H, TRUE)] <- h
    H <- H + t(H)
    diag(H) <- diag(H) / 2
    -H
}
