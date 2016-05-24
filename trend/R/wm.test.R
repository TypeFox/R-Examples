wm.test <- function(x) {
    dname <- deparse(substitute(x))
    stopifnot(is.numeric(x))
    x <- na.omit(x)
    n <- length(x)
    if(n < 2) {
        stop("sample size must be greater than 1")
    }
    x.d <- diff(x)
    x.d <- x.d[x.d != 0]
    x.s <- sign(x.d)
    h <- 0
    nn <- length(x.s)
    if (nn > 2) {
        for (i in 2:nn) {
            if (x.s[i] != x.s[i-1]) h <- h + 1
        }
    }
    h <- h - 1
    h <- max(c(0,h))
    if (n > 30) {
        z <- abs(h - (2 * n - 7) / 3) / sqrt((16 * n - 29)/ 90)
    } else {
        z <- (abs(h - (2 * n - 7) / 3) - 0.5) /
            sqrt((16 * n - 29)/ 90)
    }
    names(z) <- "Wallis and Moore z-value"
    if (z >= 0) {
        p.value <- 2 * pnorm(z, lower.tail=FALSE)
    } else {
        p.value <- 2 * pnorm(z, lower.tail=TRUE)
    }
    out <- list(statistic = z, p.value = p.value,
                alternative = "The series is significantly different from randomness",
                data.name = dname,
                method = "Wallis and Moore phase-frequency test")
    class(out) <- "htest"
    out
}
