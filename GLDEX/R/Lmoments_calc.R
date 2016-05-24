`Lmoments_calc` <-
function (data, rmax = 4) 
{
    data <- as.matrix(data)
    n <- dim(data)[1]
    p <- dim(data)[2]
    x <- array(, c(p, n))
    L <- array(, c(p, rmax))
    for (i in 1:p) {
        x[i, ] <- sort(data[, i])
    }
    if (rmax == 1) 
        return(rowMeans(x))
    bcoef <- array(, c(rmax, n))
    bcoefm <- array(, c(rmax, p, n))
    b <- array(, c(p, rmax))
    bcoef[1, ] <- seq(0, 1, by = (1/(n - 1)))
    bcoefm[1, , ] <- t(array(rep(bcoef[1, ], p), c(n, p)))
    b[, 1] <- rowMeans(x)
    b[, 2] <- rowMeans(bcoefm[1, , ] * x)
    L[, 1] = b[, 1]
    if (rmax > 2) {
        for (r in 2:(rmax - 1)) {
            rr <- r + 1
            bcoef[r, ] <- bcoef[r - 1, ] * seq((-(r - 1)/(n - 
                r)), 1, by = (1/(n - r)))
            bcoefm[r, , ] <- t(array(rep(bcoef[r, ], p), c(n, 
                p)))
            b[, rr] <- rowMeans(bcoefm[r, , ] * x)
        }
    }
    for (r in 1:(rmax - 1)) {
        L[, r + 1] <- 0
        for (k in 0:r) {
            kk <- k + 1
            L[, r + 1] <- L[, r + 1] + (-1)^(r - k) * gamma(r + 
                k + 1)/(gamma(k + 1)^2)/gamma(r - k + 1) * b[, 
                kk]
        }
    }
    return(L)
}

