ppccNormSinglyCensored <-
function (x, censored, censoring.side) 
{
    n.cen <- sum(censored)
    x.no.cen <- sort(x[!censored])
    N <- length(x)
    n <- N - n.cen
    m.tilda <- qnorm(ppoints(N, a = 3/8))
    ss.m <- sum(m.tilda^2)
    c.vec <- m.tilda/sqrt(ss.m)
    if (censoring.side == "right") 
        index <- 1:n
    else index <- (n.cen + 1):N
    cor(c.vec[index], x.no.cen)
}
