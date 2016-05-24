d3 <- function(n)
{
    d2 <- d2(n)
    e <- vector()
    for(i in 1:length(n))
    {
        int <- integrate(function(w) { w * (1 - ptukey(w, n[i], Inf)) }, 0, Inf)
        e <- append(e, sqrt(2 * int[[1]] - (d2[1]) ^ 2))
    }
    return(e)
}