Mv <-
function (u, alpha, beta, gama, low, upp, N) 
{
    x <- -log(seq(from = exp(-low), to = exp(-upp), length = N))
    f <- alpha/gamma(1 - gama) * x^(-(1 + gama)) * exp(-(u + 
        beta) * x)
    dx <- diff(x)
    h <- (f[-1] + f[-N])/2
    Mv <- rep(0, N)
    for (i in seq(N - 1, 1)) Mv[i] <- Mv[i + 1] + dx[i] * h[i]
    return(list(v = x, Mv = Mv))
}
