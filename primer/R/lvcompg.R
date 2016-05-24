`lvcompg` <-
function (t, n, parms) 
{
    r <- parms[[1]]
    a <- parms[[2]]
    dns.dt <- r * n * (1 - (a %*% n))
    return(list(c(dns.dt)))
}
