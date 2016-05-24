"sf.test" <-
function (x) 
{
    DNAME <- deparse(substitute(x))
    x <- sort(x[complete.cases(x)])
    n <- length(x)
    if ((n < 5 || n > 5000)) 
        stop("sample size must be between 5 and 5000")
    y <- qnorm(ppoints(n, a = 3/8))
    W <- cor(x, y)^2
    u <- log(n)
    v <- log(u)
    mu <- -1.2725 + 1.0521 * (v - u)
    sig <- 1.0308 - 0.26758 * (v + 2/u)
    z <- (log(1 - W) - mu)/sig
    pval <- pnorm(z, lower.tail = FALSE)
    RVAL <- list(statistic = c(W = W), p.value = pval, method = "Shapiro-Francia normality test", 
        data.name = DNAME)
    class(RVAL) <- "htest"
    return(RVAL)
}
