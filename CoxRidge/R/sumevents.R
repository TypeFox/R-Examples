sumevents <-
function (i, X, Ftime, theta, n, p, eventtimes)
{
    xj <- matrix(X[i:n, ], nrow = (n - i + 1))
    Fi <- Ftime[eventtimes == i, ]
    rr <- exp(xj %*% theta %*% Fi)
    hi <- 1/sum(rr)
    xmean <- hi * t(xj) %*% rr
    xdif <- xj - matrix(rep(xmean, n - i + 1), ncol = p, byrow = TRUE)
    XX <- t(xdif) %*% (xdif * matrix(rep(rr, p), ncol = p, byrow = FALSE))
    FF <- Fi %*% t(Fi)
    list(hi, -xmean %*% t(Fi), hi * kronecker(FF, XX))
}
