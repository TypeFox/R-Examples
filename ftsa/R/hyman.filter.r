hyman.filter = function(z)
{
    n <- length(z$x)
    ss <- (z$y[2:n] - z$y[1:(n - 1)])/(z$x[2:n] - z$x[1:(n - 1)])
    S0 <- c(ss[1], ss)
    S1 <- c(ss, ss[n - 1])
    sig <- z$b
    ind <- (S0 * S1 > 0)
    sig[ind] <- S1[ind]
    ind <- (sig >= 0)
    if (sum(ind))
        z$b[ind] <- pmin(pmax(0, z$b[ind]), 3 * pmin(abs(S0[ind]),
            abs(S1[ind])))
    ind <- !ind
    if (sum(ind))
        z$b[ind] <- pmax(pmin(0, z$b[ind]), -3 * pmin(abs(S0[ind]),
            abs(S1[ind])))
    z
}
