fuzzy.ranksum.ci <- function(x, y,
    alternative = c("two.sided", "less", "greater"),
    tol = sqrt(.Machine$double.eps), conf.level = 0.95)
{
    alternative <- match.arg(alternative)

    if (! is.numeric(x))
        stop("'x' must be numeric")
    if (! all(is.finite(x)))
        stop("'x' must be all finite")

    if (! is.numeric(y))
        stop("'y' must be numeric")
    if (! all(is.finite(y)))
        stop("'y' must be all finite")

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

    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

    nx <- length(x)
    ny <- length(y)
    N <- nx * ny
    z <- sort(as.numeric(outer(x, y, "-")))

    if (conf.level == 0 || conf.level == 1) {
        whyknots <- c(-Inf, Inf)
        kvalues <- rep(NA, 2)
        ivalues <- conf.level
    } else {

        if (alternative == "two.sided") {
            m <- qwilcox(alpha / 2, nx, ny)
            while (alpha <= 2 * pwilcox(m - 1, nx, ny))
                m <- m - 1
            while (alpha > 2 * pwilcox(m, nx, ny))
                m <- m + 1
            if (m > N / 2)
                m <- floor(N / 2)
            if (m > 0)
                whyknots <- c(z[m], z[m + 1], z[N - m], z[N - m + 1])
            else
                whyknots <- c(-Inf, z[m + 1], z[N - m], Inf)

            gdenom <- dwilcox(m, nx, ny)
            if (m < N - m)
                gdenom <- 2 * gdenom
            if (m + 1 < N - m) {
                gnumer <- 2 * pwilcox(m, nx, ny) - alpha
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
                    kvalues[i] <- 1 - fuzzy.ranksum.test(x, y,
                        alternative = alternative,
                        mu = whyknots[i], tol = 0, alpha = alpha)$reject
                }
            }
        } else {
            ##### alternative != "two.sided" #####

            m <- qwilcox(alpha, nx, ny)
            while (alpha <= pwilcox(m - 1, nx, ny))
                m <- m - 1
            while (alpha > pwilcox(m, nx, ny))
                m <- m + 1
            if (m > N)
                m <- N

            gdenom <- dwilcox(m, nx, ny)
            gnumer <- pwilcox(m, nx, ny) - alpha
            g <- gnumer / gdenom

            if (alternative == "less") {
                if (m > 0)
                    whyknots <- c(-Inf, z[N - m], z[N - m + 1])
                else
                    whyknots <- c(-Inf, z[N - m], Inf)
                ivalues <- c(1, g)
            } else {
                if (m > 0)
                    whyknots <- c(z[m], z[m + 1], Inf)
                else
                    whyknots <- c(-Inf, z[m + 1], Inf)
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
                    kvalues[i] <- 1 - fuzzy.ranksum.test(x, y,
                        alternative = alternative,
                        mu = whyknots[i], tol = 0, alpha = alpha)$reject
                }
            }
        }
    }

    method <- "Mann-Whitney-Wilcoxon rank sum test"
    foo <- list(knots = whyknots, knot.values = kvalues,
        interval.values = ivalues, alternative = alternative,
        method = method, data.name = dname, conf.level = conf.level,
        tol = tol)
    return(structure(foo, class = "fuzzyrankci"))
}

