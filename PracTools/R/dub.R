dub <- function (c1, c2, Ctot, Nh, Sh, Yh.bar)
{
    Wh <- Nh/sum(Nh)
    Ybar <- sum(Wh * Yh.bar)
    V1 <- sum(Wh * (Yh.bar - Ybar)^2)
    V2 <- sum(Wh * Sh)^2
    neyman <- Wh * Sh/sum(Wh * Sh)
    K <- (V2/V1)/(c2/c1)
    n1 <- Ctot/(c1 + c2 * sqrt(K))
    n2 <- n1 * sqrt(K)
    cost.chk <- c1 * n1 + c2 * n2
    ney.alloc <- n2 * neyman
    Vopt <- V1/n1 + V2/n2
    nsrs <- Ctot/c2
    S2 <- sum(Wh * Sh^2) + V1
    Vsrs <- S2/nsrs
    output <- structure(list(V1 = V1, V2 = V2, n1 = n1, n2 = n2, `n2/n1` = n2/n1,
        ney.alloc = ney.alloc, Vopt = Vopt, nsrs = nsrs, Vsrs = Vsrs,
        Vratio = round(Vopt/Vsrs, 2), Ctot = Ctot, cost.chk = cost.chk),
        class = "power.htest")
    output
}
