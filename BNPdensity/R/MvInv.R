MvInv <-
function (eps, u, alpha, beta, gama, N) 
{
    x <- -log(seq(from = exp(-1e-05), to = exp(-10), length = N))
    f <- alpha/gamma(1 - gama) * x^(-(1 + gama)) * exp(-(u + 
        beta) * x)
    dx <- diff(x)
    h <- (f[-1] + f[-N])/2
    Mv <- rep(0, N)
    for (i in seq(N - 1, 1)) Mv[i] <- Mv[i + 1] + dx[i] * h[i]
    err <- 1
    w <- 0
    v <- NULL
    while (err > eps) {
        w <- w + rgamma(1, 1, 1)
        v <- c(v, x[which.min(Mv > w)])
        err <- min(v)/sum(v)
    }
    return(v)
}
