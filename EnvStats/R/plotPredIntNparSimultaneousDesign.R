plotPredIntNparSimultaneousDesign <-
function (x.var = "n", y.var = "conf.level", range.x.var = NULL, 
    n = max(25, lpl.rank + n.plus.one.minus.upl.rank + 1), n.median = 1, 
    k = 1, m = ifelse(x.var == "k", ceiling(max.x), 1), r = 2, 
    rule = "k.of.m", conf.level = 0.95, pi.type = "upper", lpl.rank = ifelse(pi.type == 
        "upper", 0, 1), n.plus.one.minus.upl.rank = ifelse(pi.type == 
        "lower", 0, 1), n.max = 5000, maxiter = 1000, integrate.args.list = NULL, 
    plot.it = TRUE, add = FALSE, n.points = 100, plot.col = "black", 
    plot.lwd = 3 * par("cex"), plot.lty = 1, digits = .Options$digits, 
    cex.main = par("cex"), ..., main = NULL, xlab = NULL, ylab = NULL, 
    type = "l") 
{
    x.var <- match.arg(x.var, c("n", "conf.level", "k", "m", 
        "r"))
    y.var <- match.arg(y.var, c("conf.level", "n"))
    rule <- match.arg(rule, c("k.of.m", "CA", "Modified.CA"))
    if (x.var == y.var) 
        stop("'x.var' and 'y.var' cannot denote the same quantity")
    if (x.var == "k" && y.var == "n") 
        stop(paste("The combination of x.var=\"k\" and", "y.var=\"n\" is currently not allowed"))
    if (rule != "k.of.m" && x.var == "k") 
        stop(paste("When rule=\"", rule, "\" you cannot set x.var=\"k\"", 
            sep = ""))
    if (rule == "Modified.CA" && x.var == "m") 
        stop("When rule=\"Modified.CA\" you cannot set x.var=\"m\"")
    if (missing(range.x.var)) {
        range.x.var <- switch(x.var, n = c(ifelse(pi.type == 
            "lower", lpl.rank + 1, n.plus.one.minus.upl.rank + 
            1), 50), conf.level = c(0.5, 0.99), k = c(1, 20), 
            m = c(1, 20), r = c(1, 20))
    }
    else {
        if (is.null(range.x.var) || !all(is.finite(range.x.var)) || 
            !is.vector(range.x.var, mode = "numeric") || length(range.x.var) != 
            2) 
            stop(paste("'range.x.var' must be a numeric vector of length 2", 
                "with no missing (NA), infinite(-Inf, Inf), or undefined(NaN) values"))
    }
    min.x <- range.x.var[1]
    max.x <- range.x.var[2]
    if (min.x >= max.x) 
        stop("The second element of 'range.x.var' must be larger than the first")
    if (!is.vector(n.points, mode = "numeric") || length(n.points) != 
        1 || n.points != trunc(n.points) || n.points < 2) 
        stop("'n.points' must be an integer larger than 1")
    switch(x.var, n = {
        if (min.x < 2 || min.x != trunc(min.x)) stop(paste("When x.var=\"n\" the first element of", 
            "range.x.var must be an integer greater than 1"))
        if (max.x != trunc(max.x)) stop(paste("When x.var=\"n\" the second element of", 
            "range.x.var must be an integer greater than", "the first element of range.x.var"))
    }, conf.level = {
        if (min.x < .Machine$double.eps || min.x > 1 - .Machine$double.eps) stop(paste("When x.var=\"conf.level\" the first element of range.x.var", 
            "must be a positive number between 0 and 1"))
        if (max.x > 1 - .Machine$double.eps) stop(paste("When x.var=\"conf.level\" the second element of", 
            "range.x.var must be an positive number between 0 and 1", 
            "and greater than the first element of range.x.var"))
    }, k = {
        if (min.x < 1 || min.x != trunc(min.x)) stop(paste("When x.var=\"k\" the first element of", 
            "range.x.var must be an integer greater than 0"))
        if (max.x != trunc(max.x)) stop(paste("When x.var=\"k\" the second element of", 
            "range.x.var must be an integer greater than", "the first element of range.x.var"))
    }, m = {
        if (min.x < 1 || min.x != trunc(min.x)) stop(paste("When x.var=\"m\" the first element of", 
            "range.x.var must be an integer greater than 0"))
        if (max.x != trunc(max.x)) stop(paste("When x.var=\"m\" the second element of", 
            "range.x.var must be an integer greater than", "the first element of range.x.var"))
    }, r = {
        if (min.x < 1 || min.x != trunc(min.x)) stop(paste("When x.var=\"r\" the first element of", 
            "range.x.var must be an integer greater than 0"))
        if (max.x != trunc(max.x)) stop(paste("When x.var=\"r\" the second element of", 
            "range.x.var must be an integer greater than", "the first element of range.x.var"))
    })
    if (!is.vector(n.max, mode = "numeric") || length(n.max) != 
        1 || !is.finite(n.max) || n.max != trunc(n.max) || n.max < 
        2) 
        stop("'n.max' must be a positive integer greater than 1")
    if (!is.vector(maxiter, mode = "numeric") || length(maxiter) != 
        1 || !is.finite(maxiter) || maxiter != trunc(maxiter) || 
        maxiter < 2) 
        stop("'maxiter' must be a positive integer greater than 1")
    pi.type <- match.arg(pi.type, c("upper", "lower"))
    if (pi.type == "upper") 
        lpl.rank <- 0
    else n.plus.one.minus.upl.rank <- 0
    if (length(lpl.rank) != 1 || !is.vector(lpl.rank, mode = "numeric") || 
        !is.finite(lpl.rank) || lpl.rank != trunc(lpl.rank) || 
        lpl.rank < 0 || lpl.rank >= n.max) 
        stop("'lpl.rank' must be a non-negative integer less than 'n.max'")
    if (pi.type == "lower" & lpl.rank < 1) 
        stop("When pi.type='lower', 'lpl.rank' must be a positive integer")
    if (length(n.plus.one.minus.upl.rank) != 1 || !is.vector(n.plus.one.minus.upl.rank, 
        mode = "numeric") || !is.finite(n.plus.one.minus.upl.rank) || 
        n.plus.one.minus.upl.rank != trunc(n.plus.one.minus.upl.rank) || 
        n.plus.one.minus.upl.rank < 0 || n.plus.one.minus.upl.rank >= 
        n.max) 
        stop("'n.plus.one.minus.upl.rank' must be a non-negative integer less than 'n.max'")
    if (pi.type == "upper" & n.plus.one.minus.upl.rank < 1) 
        stop(paste("When pi.type='two.sided' or pi.type='upper',", 
            "'n.plus.one.minus.upl.rank' must be a positive integer"))
    if (x.var != "n" && y.var != "n") {
        if (is.null(n) || !is.finite(n) || !is.vector(n, mode = "numeric") || 
            length(n) != 1 || n < 2 || n != trunc(n)) 
            stop("'n' must be an integer greater than 1")
        if (pi.type == "lower" && lpl.rank >= n) 
            stop("When pi.type='lower', 'lpl.rank' must be less than 'n'")
        if (pi.type == "upper" && n.plus.one.minus.upl.rank >= 
            n) 
            stop("When pi.type='lower', 'n.plus.one.minus.upl.rank' must be less than 'n'")
    }
    if (is.null(n.median) || !is.finite(n.median) || !is.vector(n.median, 
        mode = "numeric") || length(n.median) != 1 || n.median < 
        1 || n.median != trunc(n.median) || !is.odd(n.median)) 
        stop("'n.median' must be a positive odd integer")
    if (x.var != "conf.level" && y.var != "conf.level") {
        if (is.null(conf.level) || !is.finite(conf.level) || 
            !is.vector(conf.level, mode = "numeric") || length(conf.level) != 
            1 || conf.level <= .Machine$double.eps || conf.level >= 
            1 - .Machine$double.eps) {
            stop("'conf.level' must be a scalar between 0 and 1")
        }
    }
    if (x.var != "m") {
        if (is.null(m) || !is.finite(m) || !is.vector(m, mode = "numeric") || 
            length(m) != 1 || m < 1 || m != trunc(m)) 
            stop("'m' must be a positive integer")
    }
    if (x.var != "r") {
        if (is.null(r) || !is.finite(r) || !is.vector(r, mode = "numeric") || 
            length(r) != 1 || r < 1 || r != trunc(r)) 
            stop("'r' must be a positive integer")
    }
    if (x.var != "k" && rule == "k.of.m") {
        if (is.null(k) || !is.finite(k) || !is.vector(k, mode = "numeric") || 
            length(k) != 1 || k < 1 || k != trunc(k)) 
            stop("'k' must be a positive integer")
    }
    if (x.var != "k" && x.var != "m" && rule == "k.of.m" && k > 
        m) 
        stop("'k' must be less than or equal to 'm'")
    if (x.var == "k" && max.x > m) 
        stop(paste("When x.var=\"k\" 'max.x' must be", "less than or equal to 'm'"))
    if (x.var == "m" && rule == "k.of.m" && min.x < k) 
        stop(paste("When x.var=\"m\" 'min.x' must be", "greater than or equal to 'k'"))
    pi.string <- switch(pi.type, two.sided = "(Two-Sided PI)", 
        lower = "(One-Sided Lower PI)", upper = "(One-Sided Upper PI)")
    n.string <- paste("n =", n)
    n.median.string <- ifelse(n.median == 1, "", paste("n.median = ", 
        n, ", ", sep = ""))
    conf.string <- paste("Confidence Level = ", format(100 * 
        conf.level, digits = digits), "%", sep = "")
    k.string <- paste("k =", k)
    m.string <- paste("m =", m)
    r.string <- paste("r =", r)
    rule.string <- switch(rule, k.of.m = "Based on k-of-m Rule", 
        CA = "Based on California Rule", Modified.CA = "Based on Modified California Rule")
    if (x.var != "n" & y.var != "n") {
        upl.rank <- n + 1 - n.plus.one.minus.upl.rank
        rank.string <- switch(pi.type, lower = paste("Rank(LPL) =", 
            lpl.rank), upper = paste("Rank(UPL) =", upl.rank))
    }
    else {
        rank.string <- switch(pi.type, lower = paste("Rank(LPL) =", 
            lpl.rank), upper = paste("Rank(UPL) =", ifelse(n.plus.one.minus.upl.rank == 
            1, "n", paste("n -", n.plus.one.minus.upl.rank - 
            1))))
    }
    if (plot.it) 
        gen.gp.list <- checkGraphicsPars(...)$gen.gp.list
    combo <- paste(c(x.var, y.var), collapse = " & ")
    switch(combo, `n & conf.level` = {
        x <- seq(min.x, max.x, by = ceiling((max.x - min.x + 
            1)/n.points))
        if (is.null(xlab)) xlab <- "Sample Size (n)"
        y <- predIntNparSimultaneousConfLevel(n = x, n.median = n.median, 
            k = k, m = m, r = r, rule = rule, lpl.rank = lpl.rank, 
            n.plus.one.minus.upl.rank = n.plus.one.minus.upl.rank, 
            pi.type = pi.type, integrate.args.list = integrate.args.list)
        if (is.null(ylab)) ylab <- "Confidence Level"
        line1 <- paste("Confidence Level vs. Sample Size for", 
            "Nonparametric Simultaneous Prediction Interval")
        line2 <- switch(rule, k.of.m = paste(rule.string, " with ", 
            n.median.string, k.string, ", ", m.string, ", ", 
            r.string, " ", pi.string, sep = ""), CA = paste(rule.string, 
            " with ", n.median.string, m.string, ", ", r.string, 
            " ", pi.string, sep = ""), Modified.CA = paste(rule.string, 
            " with ", n.median.string, r.string, " ", pi.string, 
            sep = ""))
        line3 <- rank.string
    }, `conf.level & n` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(xlab)) xlab <- "Confidence Level"
        y <- predIntNparSimultaneousN(n.median = n.median, k = k, 
            m = m, r = r, rule = rule, lpl.rank = lpl.rank, n.plus.one.minus.upl.rank = n.plus.one.minus.upl.rank, 
            pi.type = pi.type, conf.level = x, n.max = n.max, 
            integrate.args.list = integrate.args.list, maxiter = maxiter)
        if (is.null(ylab)) ylab <- "Sample Size (n)"
        line1 <- paste("Sample Size vs. Confidence Level for", 
            "Nonparametric Simultaneous Prediction Interval")
        line2 <- switch(rule, k.of.m = paste(rule.string, " with ", 
            n.median.string, k.string, ", ", m.string, ", ", 
            r.string, " ", pi.string, sep = ""), CA = paste(rule.string, 
            " with ", n.median.string, m.string, ", ", r.string, 
            " ", pi.string, sep = ""), Modified.CA = paste(rule.string, 
            " with ", n.median.string, r.string, " ", pi.string, 
            sep = ""))
        line3 <- rank.string
    }, `k & conf.level` = {
        x <- seq(min.x, max.x, by = ceiling((max.x - min.x + 
            1)/n.points))
        if (is.null(xlab)) xlab <- "Min # Future Obs PI Should Contain (k)"
        y <- predIntNparSimultaneousConfLevel(n = n, n.median = n.median, 
            k = x, m = m, r = r, rule = rule, lpl.rank = lpl.rank, 
            n.plus.one.minus.upl.rank = n.plus.one.minus.upl.rank, 
            pi.type = pi.type, integrate.args.list = integrate.args.list)
        if (is.null(ylab)) ylab <- "Confidence Level"
        line1 <- paste("Confidence Level vs. Min # Future Obs for", 
            "Nonparametric Simultaneous Prediction Interval")
        line2 <- paste(rule.string, " with ", n.string, ", ", 
            n.median.string, , m.string, ", ", r.string, " ", 
            pi.string, sep = "")
        line3 <- rank.string
    }, `m & conf.level` = {
        x <- seq(min.x, max.x, by = ceiling((max.x - min.x + 
            1)/n.points))
        if (is.null(xlab)) xlab <- "# Future Obs (m)"
        y <- predIntNparSimultaneousConfLevel(n = n, n.median = n.median, 
            k = k, m = x, r = r, rule = rule, lpl.rank = lpl.rank, 
            n.plus.one.minus.upl.rank = n.plus.one.minus.upl.rank, 
            pi.type = pi.type, integrate.args.list = integrate.args.list)
        if (is.null(ylab)) ylab <- "Confidence Level"
        line1 <- paste("Confidence Level vs. # Future Obs for", 
            "Nonparametric Simultaneous Prediction Interval")
        line2 <- paste(rule.string, " with ", n.string, ", ", 
            n.median.string, , k.string, ", ", r.string, " ", 
            pi.string, sep = "")
        line3 <- rank.string
    }, `m & n` = {
        x <- seq(min.x, max.x, by = ceiling((max.x - min.x + 
            1)/n.points))
        if (is.null(xlab)) xlab <- "# Future Obs (m)"
        y <- predIntNparSimultaneousN(n.median = n.median, k = k, 
            m = x, r = r, rule = rule, lpl.rank = lpl.rank, n.plus.one.minus.upl.rank = n.plus.one.minus.upl.rank, 
            pi.type = pi.type, conf.level = conf.level, n.max = n.max, 
            integrate.args.list = integrate.args.list, maxiter = maxiter)
        if (is.null(ylab)) ylab <- "Sample Size (n)"
        line1 <- paste("Sample Size vs. # Future Obs for", "Nonparametric Simultaneous Prediction Interval")
        line2 <- paste(rule.string, " with ", n.median.string, 
            k.string, ", ", r.string, ", ", conf.string, " ", 
            pi.string, sep = "")
        line3 <- rank.string
    }, `r & conf.level` = {
        x <- seq(min.x, max.x, by = ceiling((max.x - min.x + 
            1)/n.points))
        if (is.null(xlab)) xlab <- "# Future Sampling Occasions (r)"
        y <- predIntNparSimultaneousConfLevel(n = n, n.median = n.median, 
            k = k, m = m, r = x, rule = rule, lpl.rank = lpl.rank, 
            n.plus.one.minus.upl.rank = n.plus.one.minus.upl.rank, 
            pi.type = pi.type, integrate.args.list = integrate.args.list)
        if (is.null(ylab)) ylab <- "Confidence Level"
        line1 <- paste("Confidence Level vs. # Future Sampling Occasions for", 
            "Nonparametric Simultaneous Prediction Interval")
        line2 <- paste(rule.string, " with ", n.string, ", ", 
            n.median.string, , k.string, ", ", m.string, " ", 
            pi.string, sep = "")
        line3 <- rank.string
    }, `r & n` = {
        x <- seq(min.x, max.x, by = ceiling((max.x - min.x + 
            1)/n.points))
        if (is.null(xlab)) xlab <- "# Future Sampling Occasions (r)"
        y <- predIntNparSimultaneousN(n.median = n.median, k = k, 
            m = m, r = x, rule = rule, lpl.rank = lpl.rank, n.plus.one.minus.upl.rank = n.plus.one.minus.upl.rank, 
            pi.type = pi.type, conf.level = conf.level, n.max = n.max, 
            integrate.args.list = integrate.args.list, maxiter = maxiter)
        if (is.null(ylab)) ylab <- "Sample Size (n)"
        line1 <- paste("Sample Size vs. # Future Sampling Occasions for", 
            "Nonparametric Simultaneous Prediction Interval")
        line2 <- paste(rule.string, " with ", n.median.string, 
            k.string, ", ", m.string, ", ", conf.string, " ", 
            pi.string, sep = "")
        line3 <- rank.string
    })
    if (plot.it) {
        if (!add) {
            plot(x, y, type = "n", main = "", sub = "", ..., 
                xlab = xlab, ylab = ylab)
            if (is.null(main)) {
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
    names(ret.list) <- c(x.var, y.var)
    invisible(ret.list)
}
