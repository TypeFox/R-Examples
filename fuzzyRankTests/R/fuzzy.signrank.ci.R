fuzzy.signrank.ci <- function(x,
    alternative = c("two.sided", "less", "greater"),
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

    nx <- length(x)
    N <- nx * (nx + 1) / 2
    w <- outer(x, x, "+") / 2
    w <- sort(w[lower.tri(w, diag = TRUE)])

    if (conf.level == 0 || conf.level == 1) {
        whyknots <- c(-Inf, Inf)
        kvalues <- rep(NA, 2)
        ivalues <- conf.level
    } else {

        if (alternative == "two.sided") {
            m <- qsignrank(alpha / 2, nx)
            while (alpha <= 2 * psignrank(m - 1, nx))
                m <- m - 1
            while (alpha > 2 * psignrank(m, nx))
                m <- m + 1
            if (m > N / 2)
                m <- floor(N / 2)
            if (m > 0)
                whyknots <- c(w[m], w[m + 1], w[N - m], w[N - m + 1])
            else
                whyknots <- c(-Inf, w[m + 1], w[N - m], Inf)

            gdenom <- dsignrank(m, nx)
            if (m < N - m)
                gdenom <- 2 * gdenom
            if (m + 1 < N - m) {
                gnumer <- 2 * psignrank(m, nx) - alpha
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
                    kvalues[i] <- 1 - fuzzy.signrank.test(x,
                        alternative = alternative,
                        mu = whyknots[i], tol = 0, alpha = alpha)$reject
                }
            }
        } else {
            ##### alternative != "two.sided" #####

            m <- qsignrank(alpha, nx)
            while (alpha <= psignrank(m - 1, nx))
                m <- m - 1
            while (alpha > psignrank(m, nx))
                m <- m + 1
            if (m > N)
                m <- N

            gdenom <- dsignrank(m, nx)
            gnumer <- psignrank(m, nx) - alpha
            g <- gnumer / gdenom

            if (alternative == "less") {
                if (m > 0)
                    whyknots <- c(-Inf, w[N - m], w[N - m + 1])
                else
                    whyknots <- c(-Inf, w[N - m], Inf)
                ivalues <- c(1, g)
            } else {
                if (m > 0)
                    whyknots <- c(w[m], w[m + 1], Inf)
                else
                    whyknots <- c(-Inf, w[m + 1], Inf)
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
                    kvalues[i] <- 1 - fuzzy.signrank.test(x,
                        alternative = alternative,
                        mu = whyknots[i], tol = 0, alpha = alpha)$reject
                }
            }
        }
    }

    method <- "Wilcoxon signed rank test"
    foo <- list(knots = whyknots, knot.values = kvalues,
        interval.values = ivalues, alternative = alternative,
        method = method, data.name = dname, conf.level = conf.level,
        tol = tol)
    return(structure(foo, class = "fuzzyrankci"))
}

