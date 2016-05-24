chenTTest.sub <-
function (mu, muhat, sdhat, skewhat, n, alternative) 
{
    a <- skewhat/(6 * sqrt(n))
    t.stat <- (muhat - mu)/(sdhat/sqrt(n))
    t.stat.chen <- t.stat + a * (1 + 2 * t.stat^2) + 4 * a^2 * 
        (t.stat + 2 * t.stat^3)
    names(t.stat.chen) <- "t"
    df <- n - 1
    names(df) <- "df"
    switch(alternative, greater = {
        if (skewhat <= 0) {
            warning(paste("The sample skew is less than or equal to 0.\n ", 
                "Chen's test is not appropriate for a\n ", "\"greater\" alternative.\n"))
            p.z <- NA
            p.t <- NA
            p.avg <- NA
        } else {
            p.z <- 1 - pnorm(t.stat.chen)
            p.t <- 1 - pt(t.stat.chen, df = df)
            p.avg <- mean(c(p.z, p.t))
        }
    }, less = {
        if (skewhat >= 0) {
            warning(paste("The sample skew is greater than or equal to 0.\n ", 
                "Chen's test is not appropriate for a\n ", "\"less\" alternative.\n"))
            p.z <- NA
            p.t <- NA
            p.avg <- NA
        } else {
            p.z <- pnorm(t.stat.chen)
            p.t <- pt(t.stat.chen, df = df)
            p.avg <- mean(c(p.z, p.t))
        }
    })
    p.value <- c(p.z, p.t, p.avg)
    names(p.value) <- c("z", "t", "Avg. of z and t")
    list(statistic = t.stat.chen, parameters = df, p.value = p.value)
}
