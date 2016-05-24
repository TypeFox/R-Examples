fci.binom <- function(x, n, alpha = 0.05, p = seq(0, 1, length = 10001),
    flat = 1 / 4) {

    if (! is.numeric(x)) stop("x not numeric")
    if (! is.numeric(n)) stop("n not numeric")
    if (! is.numeric(alpha)) stop("alpha not numeric")
    if (! is.numeric(p)) stop("p not numeric")

    if (length(x) != 1) stop("x not scalar")
    if (length(n) != 1) stop("n not scalar")
    if (length(alpha) != 1) stop("alpha not scalar")

    if (as.integer(x) != x) stop("x not integer")
    if (as.integer(n) != n) stop("n not integer")

    if (! (n > 0)) stop("n not positive")
    if (! (0 <= x & x <= n)) stop("x not in 0, ..., n")
    if (! (0 < alpha & alpha < 1)) stop("alpha not in (0, 1)")
    if (! all(0 <= p & p <= 1)) stop("p not in [0, 1]")

    phi <- umpu.binom(x, n, p, alpha)
    support <- range(p[phi < 1])

    cat(100 * (1 - alpha), "percent fuzzy confidence interval\n")

    if (x == 0 || x == n) {

        cat("core is empty\n")
        if (x == 0) {
            cat("support is [", support[1], ", ", support[2], ")\n", sep = "")
            xlim <- c(0, support[2] * (1 + flat))
        } else {
            cat("support is (", support[1], ", ", support[2], "]\n", sep = "")
            xlim <- c(support[1] - (1 - support[1]) * flat, 1)
        }

        plot(p, 1 - phi, xlim = xlim, ylim = c(0, 1),
            ylab = expression(1 - phi(x, alpha, p)), type = "l")

    } else {

        core <- range(p[phi == 0])
        cat("core is [", core[1], ", ", core[2], "]\n", sep = "")
        cat("support is (", support[1], ", ", support[2], ")\n", sep = "")

        lim.low <- c(support[1], core[1])
        lim.hig <- c(core[2], support[2])

        extra <- flat * max(diff(lim.low), diff(lim.hig))
        lim.low <- range(lim.low - extra, lim.low + extra)
        lim.hig <- range(lim.hig - extra, lim.hig + extra)
        lim.low <- pmax(0, lim.low)
        lim.hig <- pmax(0, lim.hig)
        lim.low <- pmin(1, lim.low)
        lim.hig <- pmin(1, lim.hig)

        if (lim.low[2] < lim.hig[1]) {
            oldpar <- par(mfrow = c(1, 2))
            plot(p, 1 - phi, xlim = lim.low, ylim = c(0, 1),
                ylab = expression(1 - phi(x, alpha, p)), type = "l")
            plot(p, 1 - phi, xlim = lim.hig, ylim = c(0, 1),
                ylab = expression(1 - phi(x, alpha, p)), type = "l")
            par(mfrow = oldpar)
        } else {
            plot(p, 1 - phi, xlim = c(lim.low[1], lim.hig[2]), ylim = c(0, 1),
                ylab = expression(1 - phi(x, alpha, p)), type = "l")
        }
    }
}

