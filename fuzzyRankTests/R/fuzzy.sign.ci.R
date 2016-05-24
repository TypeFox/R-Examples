fuzzy.sign.ci <- function(x, alternative = c("two.sided", "less", "greater"),
    tol = sqrt(.Machine$double.eps), conf.level = 0.95)
{
    alternative <- match.arg(alternative)

    if (! is.numeric(x))
        stop("'x' must be numeric")
    if (! all(is.finite(x)))
        stop("'x' must be all finite")

    if (! is.numeric(tol))
        stop("'tol' must be numeric")
    if (length(tol) != 1)
        stop("'tol' must be a single number")
    if (! is.finite(tol))
        stop("'tol' must be finite")
    if (tol < 0.0)
        stop("'tol' must be nonnegative")

    if (! is.numeric(conf.level))
        stop("'conf.level' must be numeric")
    if (length(conf.level) != 1)
        stop("'conf.level' must be a single number")
    if (! (is.finite(conf.level) & 0 <= conf.level & conf.level <= 1))
        stop("'conf.level' must satisfy 0 <= conf.level <= 1")
    alpha <- 1 - conf.level

    dname <- deparse(substitute(x))

    n <- length(x)
    x <- sort(x)

    if (conf.level == 0 || conf.level == 1) {
        whyknots <- c(-Inf, Inf)
        kvalues <- rep(NA, 2)
        ivalues <- conf.level
    } else {

        if (alternative == "two.sided") {
            m <- qbinom(alpha / 2, n, 1 / 2)
            while (alpha <= 2 * pbinom(m - 1, n, 1 / 2))
                m <- m - 1
            while (alpha > 2 * pbinom(m, n, 1 / 2))
                m <- m + 1
            if (m > n / 2)
                m <- floor(n / 2)
            if (m > 0)
                whyknots <- c(x[m], x[m + 1], x[n - m], x[n - m + 1])
            else
                whyknots <- c(-Inf, x[m + 1], x[n - m], Inf)

            gdenom <- dbinom(m, n, 1 / 2)
            if (m < n - m)
                gdenom <- 2 * gdenom
            if (m + 1 < n - m) {
                gnumer <- 2 * pbinom(m, n, 1 / 2) - alpha
            } else {
                gnumer <- conf.level
            }
            g <- gnumer / gdenom

            ivalues <- c(g, 1, g)
            i <- 1
            while (i <= length(ivalues)) {
                if (whyknots[i] + tol >= whyknots[i + 1]) {
                    whyknots <- whyknots[- (i + 1)]
                    ivalues <- ivalues[- i]
                } else {
                    i <- i + 1
                }
            }

            kvalues <- rep(NA, length(whyknots))
            for (i in 1:length(whyknots)) {
                if (! is.finite(whyknots[i])) {
                    kvalues[i] <- NA
                } else {
                    kvalues[i] <- 1 - fuzzy.sign.test(x,
                        alternative = alternative,
                        mu = whyknots[i], tol = 0, alpha = alpha)$reject
                }
            }
        } else {
            ##### alternative != "two.sided" #####

            m <- qbinom(alpha, n, 1 / 2)
            while (alpha <= pbinom(m - 1, n, 1 / 2))
                m <- m - 1
            while (alpha > pbinom(m, n, 1 / 2))
                m <- m + 1
            if (m > n)
                m <- n

            gdenom <- dbinom(m, n, 1 / 2)
            gnumer <- pbinom(m, n, 1 / 2) - alpha
            g <- gnumer / gdenom

            if (alternative == "less") {
                if (m > 0)
                    whyknots <- c(-Inf, x[n - m], x[n - m + 1])
                else
                    whyknots <- c(-Inf, x[n - m], Inf)
                ivalues <- c(1, g)
            } else {
                if (m > 0)
                    whyknots <- c(x[m], x[m + 1], Inf)
                else
                    whyknots <- c(-Inf, x[m + 1], Inf)
                ivalues <- c(g, 1)
            }

            i <- 1
            while (i <= length(ivalues)) {
                if (whyknots[i] + tol >= whyknots[i + 1]) {
                    whyknots <- whyknots[- (i + 1)]
                    ivalues <- ivalues[- i]
                } else {
                    i <- i + 1
                }
            }

            kvalues <- rep(NA, length(whyknots))
            for (i in 1:length(whyknots)) {
                if (! is.finite(whyknots[i])) {
                    kvalues[i] <- NA
                } else {
                    kvalues[i] <- 1 - fuzzy.sign.test(x,
                        alternative = alternative,
                        mu = whyknots[i], tol = 0, alpha = alpha)$reject
                }
            }
        }
    }

    method <- "sign test"
    foo <- list(knots = whyknots, knot.values = kvalues,
        interval.values = ivalues, alternative = alternative,
        method = method, data.name = dname, conf.level = conf.level,
        tol = tol)
    return(structure(foo, class = "fuzzyrankci"))
}

