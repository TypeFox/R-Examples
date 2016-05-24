d2 <- function(n)
{
    d <- vector()
    for(i in 1:length(n))
    {
        int <- integrate(function(w) {1 - ptukey(w, n[i], Inf)}, 0, Inf)
        d <- append(d, int[[1]])
    }
    return(d)
}