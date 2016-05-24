HessLongPH <-
function (X, Xtime, Xtime2, ew) {
    H1 <- - XtX / sigma^2
    H2 <- matrix(0, ncx, ncx)
    index <- which(lower.tri(H2, TRUE), arr.ind = TRUE)
    nn <- nrow(index)
    h <- numeric(nn)
    for (i in 1:nn) {
        i1 <- index[i, 1]
        i2 <- index[i, 2]
        pp <- rowsum(lambda0. * (alpha^2 * Xtime2[, i1] * Xtime2[, i2]) * ew, indT)
        h[i] <- - sum((pp * p.byt.) %*% wGH)
    }
    H2[lower.tri(H2, TRUE)] <- h
    H2 <- H2 + t(H2)
    diag(H2) <- diag(H2) / 2
    -(H1 + H2)
}
