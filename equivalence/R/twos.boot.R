"twos.boot" <-
function (xy, i, Epsilon=0.5)
{
    x <- xy[i, 1]
    y <- xy[i, 2]
    n1 <- length(x)
    n2 <- length(y)
    mean <- mean(y) - mean(x)
    sp <- sqrt(((n1-1)*var(x) 
              + (n2-1)*var(y))/(n1+n2-2))

    df1 <- 1
    df2 <- n1+n2-2
    cutoff <- c(0)
    ncp <- Epsilon^2 * n1*n2/(n1+n2-2)

    findquant <- function(alpha, q, df1, df2, ncp) alpha - pf(q,
        df1, df2, ncp)
    cutoff <- sqrt(uniroot(findquant, c(0, 1e+23), alpha = 0.05,
        df1 = df1, df2 = min(df2, 4000), ncp = ncp)$root)
    tstat <- mean/(sp*sqrt(1/n1+1/n2))
    result <- (abs(tstat) < cutoff)
}

