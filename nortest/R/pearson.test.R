"pearson.test" <-
function (x, n.classes = ceiling(2 * (n^(2/5))), adjust = TRUE) 
{
    DNAME <- deparse(substitute(x))
    x <- x[complete.cases(x)]
    n <- length(x)
    if (adjust) {
        dfd <- 2
    }
    else {
        dfd <- 0
    }
    num <- floor(1 + n.classes * pnorm(x, mean(x), sd(x)))
    count <- tabulate(num, n.classes)
    prob <- rep(1/n.classes, n.classes)
    xpec <- n * prob
    h <- ((count - xpec)^2)/xpec
    P <- sum(h)
    pvalue <- pchisq(P, n.classes - dfd - 1, lower.tail = FALSE)
    RVAL <- list(statistic = c(P = P), p.value = pvalue, method = "Pearson chi-square normality test", 
        data.name = DNAME, n.classes = n.classes, df = n.classes - 
            1 - dfd)
    class(RVAL) <- "htest"
    return(RVAL)
}
