gs3 <-
function (ut, n, r, alpha, beta, gama, delta) 
{
    w <- ut
    ratio <- NaN
    while (is.nan(ratio)) {
        v <- ustar <- rgamma(1, shape = delta, rate = delta/ut)
        vw <- v/w
        vb <- v + beta
        wb <- w + beta
        A <- vw^(n - 2 * delta)
        B <- (vb/wb)^(r * gama - n)
        D <- vb^gama - wb^gama
        E <- 1/vw - vw
        ratio <- A * B * exp(-alpha/gama * D - delta * E)
    }
    p <- min(1, ratio)
    u <- ifelse(runif(1) <= p, ustar, ut)
    return(u)
}
