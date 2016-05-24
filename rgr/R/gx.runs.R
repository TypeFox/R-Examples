gx.runs <-
function (n1, n2, u) 
{
    m <- 2 * n1 * n2
    n <- n1 + n2
    Eu <- 1 + m/n
    Varu <- (m * (m - n))/(n * n * (n - 1))
    z <- (u - Eu)/sqrt(Varu)
    pz <- 100 * round((pnorm(z)), 3)
    cat("  Traverse length is", n, "sites with", n2, "anomalous sites in", 
        u, "runs\n  E(u) =", round(Eu, 2), ", Var(u) =", round(Varu, 
            2), "and z =", signif(z, 4), paste("\n  Probability that this is due to 'chance' is ", 
            pz, "%\n", sep = ""))
    invisible()
}
