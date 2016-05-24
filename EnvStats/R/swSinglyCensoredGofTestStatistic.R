swSinglyCensoredGofTestStatistic <-
function (x, censored, censoring.side) 
{
    n.cen <- sum(censored)
    x.no.cen <- sort(x[!censored])
    N <- length(x)
    n <- N - n.cen
    m.tilda <- qnorm(ppoints(N, a = 3/8))
    ss.m <- sum(m.tilda^2)
    c.vec <- m.tilda/sqrt(ss.m)
    a.tilda <- numeric(N)
    y <- 1/sqrt(N)
    a.tilda[N] <- c.vec[N] + 0.221157 * y - 0.147981 * y^2 - 
        2.07119 * y^3 + 4.434685 * y^4 - 2.706056 * y^5
    a.tilda[1] <- -a.tilda[N]
    if ((2 < N) & (N <= 5)) {
        phi <- (ss.m - 2 * m.tilda[N]^2)/(1 - 2 * a.tilda[N]^2)
        a.tilda[2:(N - 1)] <- m.tilda[2:(N - 1)]/sqrt(phi)
    }
    else {
        a.tilda[N - 1] <- c.vec[N - 1] + 0.042981 * y - 0.293762 * 
            y^2 - 1.752461 * y^3 + 5.682633 * y^4 - 3.582663 * 
            y^5
        phi <- (ss.m - 2 * m.tilda[N]^2 - 2 * m.tilda[N - 1]^2)/(1 - 
            2 * a.tilda[N]^2 - 2 * a.tilda[N - 1]^2)
        a.tilda[3:(N - 2)] <- m.tilda[3:(N - 2)]/sqrt(phi)
        a.tilda[2] <- -a.tilda[N - 1]
    }
    if (censoring.side == "right") 
        index <- 1:n
    else index <- (n.cen + 1):N
    cor(a.tilda[index], x.no.cen)^2
}
