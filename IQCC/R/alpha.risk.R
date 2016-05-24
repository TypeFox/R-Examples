alpha.risk <- function(n)
{
    D1 <- function(n)
    {
        d1 <- max(0, d2(n) - 3 * d3(n))
        return(d1)
    }
    D2 <- function(n)
    {
        D2 <- d2(n) + 3 * d3(n)
        return(D2)
    }
    risco <- function(n)
    {
        risco <- 1 - (ptukey(D2(n), n, Inf) - ptukey(D1(n), n, Inf))
        return(risco)
    }
    risk <- rep(0, length(n))
    for(i in 1:length(n))
        risk[i] <- risco(n[i])
    return(risk)
}