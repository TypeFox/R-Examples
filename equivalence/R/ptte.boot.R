"ptte.boot" <-
function (x, i, Epsilon=0.5)
{
    x <- x[i]
    mean <- mean(x)
    std <- sd(x)
    n = length(x)

    df1 <- 1
    df2 <- n - 1
    cutoff <- c(0)
    ncp <- Epsilon^2 * n

    findquant <- function(alpha, q, df1, df2, ncp) alpha - pf(q,
        df1, df2, ncp)
    cutoff <- sqrt(uniroot(findquant, c(0, 1e+23), alpha = 0.05,
        df1 = df1, df2 = min(df2, 4000), ncp = ncp)$root)
    tstat <- mean/(std/sqrt(n))
    result <- (abs(tstat) < cutoff)
}

