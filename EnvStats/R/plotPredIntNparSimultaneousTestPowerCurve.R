plotPredIntNparSimultaneousTestPowerCurve <-
function (n = 8, n.median = 1, k = 1, m = 2, r = 1, rule = "k.of.m", 
    lpl.rank = ifelse(pi.type == "upper", 0, 1), n.plus.one.minus.upl.rank = ifelse(pi.type == 
        "lower", 0, 1), pi.type = "upper", r.shifted = r, integrate.args.list = NULL, 
    method = "approx", NMC = 100, range.delta.over.sigma = c(0, 
        5), plot.it = TRUE, add = FALSE, n.points = 20, plot.col = "black", 
    plot.lwd = 3 * par("cex"), plot.lty = 1, digits = .Options$digits, 
    cex.main = par("cex"), ..., main = NULL, xlab = NULL, ylab = NULL, 
    type = "l") 
{
    rule <- match.arg(rule, c("k.of.m", "CA", "Modified.CA"))
    pi.type <- match.arg(pi.type, c("upper", "lower"))
    method <- match.arg(method, c("approx", "Monte.Carlo"))
    if (is.null(range.delta.over.sigma) || !all(is.finite(range.delta.over.sigma)) || 
        !is.vector(range.delta.over.sigma, mode = "numeric") || 
        length(range.delta.over.sigma) != 2) 
        stop(paste("'range.delta.over.sigma' must be a numeric vector of length 2", 
            "with no missing (NA), infinite(-Inf, Inf), or undefined(NaN) values"))
    min.x <- range.delta.over.sigma[1]
    max.x <- range.delta.over.sigma[2]
    if (min.x >= max.x) 
        stop("The second element of 'range.delta.over.sigma' must be larger than the first")
    if (!is.vector(n.points, mode = "numeric") || length(n.points) != 
        1 || n.points != trunc(n.points) || n.points < 2) 
        stop("'n.points' must be an integer larger than 1")
    if (is.null(n) || !is.finite(n) || !is.vector(n, mode = "numeric") || 
        length(n) != 1 || n < 2 || n != trunc(n)) 
        stop("'n' must be an integer greater than 1")
    if (is.null(n.median) || !is.finite(n.median) || !is.vector(n.median, 
        mode = "numeric") || length(n.median) != 1 || n.median < 
        1 || n.median != trunc(n.median) || !is.odd(n.median)) 
        stop("'n.median' must be a positive odd integer")
    if (rule == "k.of.m") {
        if (length(m) != 1 || !is.vector(m, mode = "numeric") || 
            !is.finite(m) || m < 1 || m != trunc(m)) 
            stop("'m' must be a positive integer")
        if (length(k) != 1 || !is.vector(k, mode = "numeric") || 
            !is.finite(k) || k < 1 || k > m || k != trunc(k)) 
            stop("'k' must be a positive integer between 1 and 'm'")
    }
    else if (rule == "CA") {
        if (length(m) != 1 || !is.vector(m, mode = "numeric") || 
            !is.finite(m) || m < 1 || m != trunc(m)) 
            stop("'m' must be a positive integer")
        if (m == 1) 
            rule <- "k.of.m"
    }
    else {
        m <- 4
    }
    if (length(r) != 1 || !is.vector(r, mode = "numeric") || 
        !is.finite(r) || r != trunc(r) || r < 1) 
        stop("'r' must be a positive integer")
    if (is.null(r.shifted) || !is.finite(r.shifted) || !is.vector(r.shifted, 
        mode = "numeric") || length(r.shifted) != 1 || r.shifted != 
        trunc(r.shifted) || r.shifted < 1 || r.shifted > r) 
        stop("'r.shifted' must be a positive integer less than or equal to 'r'")
    if (pi.type == "upper") 
        lpl.rank <- 0
    else n.plus.one.minus.upl.rank <- 0
    if (!is.vector(lpl.rank, mode = "numeric") || length(lpl.rank) != 
        1 || !is.finite(lpl.rank) || lpl.rank != trunc(lpl.rank) || 
        lpl.rank < 0 || lpl.rank >= n) 
        stop("'lpl.rank' must be a non-negative integers less than 'n'")
    if (pi.type == "lower" & lpl.rank < 1) 
        stop("When pi.type='lower', 'lpl.rank' must be a positive integer")
    if (!is.vector(n.plus.one.minus.upl.rank, mode = "numeric") || 
        length(n.plus.one.minus.upl.rank) != 1 || !is.finite(n.plus.one.minus.upl.rank) || 
        n.plus.one.minus.upl.rank != trunc(n.plus.one.minus.upl.rank) || 
        n.plus.one.minus.upl.rank < 0 || n.plus.one.minus.upl.rank >= 
        n) 
        stop("'n.plus.one.minus.upl.rank' must be a non-negative integer less than 'n'")
    if (pi.type == "upper" & n.plus.one.minus.upl.rank < 1) 
        stop("When pi.type='upper', 'n.plus.one.minus.upl.rank' must be a positive integer")
    pl.rank <- ifelse(pi.type == "upper", n + 1 - n.plus.one.minus.upl.rank, 
        lpl.rank)
    delta.string <- "(Future Mean - Background Mean) / SD"
    pi.string <- switch(pi.type, lower = "Lower One-Sided PI", 
        upper = "Upper One-Sided PI")
    n.string <- paste("n =", n)
    n.median.string <- ifelse(n.median == 1, "", paste("n.median = ", 
        n.median, ", ", sep = ""))
    pl.string <- paste("PL Rank =", pl.rank)
    if (rule == "k.of.m") 
        k.string <- paste("k =", k)
    m.string <- paste("m =", m)
    r.string <- paste("r = ", r)
    if (r > 1) 
        r.string <- paste(r.string, ", r shifted = ", r.shifted, 
            sep = "")
    rule.string <- switch(rule, k.of.m = "k-of-m", CA = "CA", 
        Modified.CA = "Modified CA")
    method.string <- switch(method, approx = "Approximating K-multiplier", 
        Monte.Carlo = paste(NMC, "Monte Carlo Simulations"))
    if (plot.it) 
        gen.gp.list <- checkGraphicsPars(...)$gen.gp.list
    x <- seq(min.x, max.x, length = n.points)
    y <- predIntNparSimultaneousTestPower(n = n, n.median = n.median, 
        k = k, m = m, r = r, rule = rule, lpl.rank = lpl.rank, 
        n.plus.one.minus.upl.rank = n.plus.one.minus.upl.rank, 
        delta.over.sigma = x, pi.type = pi.type, r.shifted = r.shifted, 
        method = method, NMC = NMC, integrate.args.list = integrate.args.list)
    if (is.null(xlab)) 
        xlab <- delta.string
    if (is.null(ylab)) 
        ylab <- "Power"
    if (plot.it) {
        if (!add) {
            plot(x, y, type = "n", main = "", sub = "", ..., 
                xlab = xlab, ylab = ylab)
            if (is.null(main)) {
                string <- switch(rule, k.of.m = paste(k.string, 
                  m.string, r.string, sep = ", "), CA = paste(m.string, 
                  r.string, sep = ", "), Modified.CA = r.string)
                line1 <- paste("Power vs. Delta/Sigma for", "Simultaneous Nonparametric Prediction Interval")
                line2 <- paste("Using ", rule.string, " Rule with ", 
                  n.string, ", ", n.median.string, pl.string, 
                  ", ", string, sep = "")
                line3 <- paste(pi.string, "Based on", method.string, 
                  "and Assuming a Normal Distribution")
                mtext(text = line1, side = 3, line = 3, cex = cex.main)
                mtext(text = line2, side = 3, line = 2, cex = cex.main)
                mtext(text = line3, side = 3, line = 1, cex = cex.main)
            }
            else {
                arg.list <- c(list(main = main), gen.gp.list, 
                  list(cex = cex.main))
                do.call("title", arg.list)
            }
            arg.list <- c(list(x = x, y = y), gen.gp.list, list(type = type, 
                col = plot.col, lwd = plot.lwd, lty = plot.lty))
            do.call("lines", arg.list)
        }
        else {
            arg.list <- c(list(x = x, y = y), gen.gp.list, list(type = type, 
                col = plot.col, lwd = plot.lwd, lty = plot.lty))
            do.call("lines", arg.list)
        }
    }
    ret.list <- list(x, y)
    names(ret.list) <- c("delta.over.sigma", "power")
    invisible(ret.list)
}
