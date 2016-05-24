swGofTestStatistic <-
function (x) 
{
    x <- sort(x)
    n <- length(x)
    m.tilda <- qnorm(ppoints(n, a = 3/8))
    ss.m <- sum(m.tilda^2)
    c.vec <- m.tilda/sqrt(ss.m)
    a.tilda <- numeric(n)
    y <- 1/sqrt(n)
    a.tilda[n] <- c.vec[n] + 0.221157 * y - 0.147981 * y^2 - 
        2.07119 * y^3 + 4.434685 * y^4 - 2.706056 * y^5
    a.tilda[1] <- -a.tilda[n]
    if (2 < n & n <= 5) {
        phi <- (ss.m - 2 * m.tilda[n]^2)/(1 - 2 * a.tilda[n]^2)
        a.tilda[2:(n - 1)] <- m.tilda[2:(n - 1)]/sqrt(phi)
    }
    else {
        a.tilda[n - 1] <- c.vec[n - 1] + 0.042981 * y - 0.293762 * 
            y^2 - 1.752461 * y^3 + 5.682633 * y^4 - 3.582663 * 
            y^5
        phi <- (ss.m - 2 * m.tilda[n]^2 - 2 * m.tilda[n - 1]^2)/(1 - 
            2 * a.tilda[n]^2 - 2 * a.tilda[n - 1]^2)
        a.tilda[3:(n - 2)] <- m.tilda[3:(n - 2)]/sqrt(phi)
        a.tilda[2] <- -a.tilda[n - 1]
    }
    (sum(a.tilda * x)^2)/((n - 1) * var(x))
}
