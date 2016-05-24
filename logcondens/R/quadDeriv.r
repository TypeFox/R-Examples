quadDeriv <- function(dx, w, eta){
    n <- length(w)
    k <- 2:n
    a <- c(0, (dx * eta)[-1])
    b <- exp(cumsum(a))
    dL2n <- 1:n * NA
    H0 <- exp(eta[1])
    dL1 <- sum(w) - H0 * sum(1/eta[2:n] * (b[2:n] - b[(2:n) - 1]))
    temp <- (dx * rev(cumsum(rev(c(0, 1/eta[3:n] * (b[3:n] - b[2:(n - 1)]), 0)))))[k]
    dL2n[k] <- rev(cumsum(rev(w[k]))) * dx[k] - H0 * (-eta[k]^(-2) * (b[k] - b[k - 1]) + 1/eta[k] * (dx[k] * b[k]) + temp)
    G <- c(dL1, dL2n[2:n])
    d22n <- 1:n * NA
    d21 <- -sum(w) + dL1
    temp <- (dx^2 * rev(cumsum(rev(c(0, 1/eta[3:n] * (b[3:n] - b[2:(n - 1)]), 0)))))[k]
    d22n[k] <- -H0 * (2 * eta[k]^(-3) * (b[k] - b[k - 1]) - 2 * eta[k]^(-2) * (dx[k] * b[k]) + dx[k]^2 * b[k]/eta[k] + temp)
    H <- c(d21, d22n[2:n])
    H[H > 0] <- -H[H > 0]
    return(-cbind(G, H))
}
