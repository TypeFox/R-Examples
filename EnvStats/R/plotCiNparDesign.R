plotCiNparDesign <-
function (x.var = "n", y.var = "conf.level", range.x.var = NULL, 
    n = 25, p = 0.5, conf.level = 0.95, ci.type = "two.sided", 
    lcl.rank = ifelse(ci.type == "upper", 0, 1), n.plus.one.minus.ucl.rank = ifelse(ci.type == 
        "lower", 0, 1), plot.it = TRUE, add = FALSE, n.points = 100, 
    plot.col = "black", plot.lwd = 3 * par("cex"), plot.lty = 1, 
    digits = .Options$digits, cex.main = par("cex"), ..., main = NULL, 
    xlab = NULL, ylab = NULL, type = "l") 
{
    x.var <- match.arg(x.var, c("n", "p", "conf.level"))
    y.var <- match.arg(y.var, c("conf.level", "n"))
    ci.type <- match.arg(ci.type, c("two.sided", "lower", "upper"))
    if (x.var == y.var) 
        stop("'x.var' and 'y.var' cannot denote the same quantity")
    if (missing(range.x.var)) {
        range.x.var <- switch(x.var, n = c(2, 50), p = c(0.5, 
            0.99), conf.level = c(0.5, 0.99))
    }
    else {
        if (is.null(range.x.var) || !all(is.finite(range.x.var)) || 
            !is.vector(range.x.var, mode = "numeric") || length(range.x.var) != 
            2 || range.x.var[1] >= range.x.var[2]) 
            stop(paste("'range.x.var' must be a numeric vector of length 2", 
                "with no missing (NA), infinite(-Inf, Inf), or undefined(NaN) values", 
                "and the first element must be less than the second element"))
    }
    min.x <- range.x.var[1]
    max.x <- range.x.var[2]
    if (!is.vector(n.points, mode = "numeric") || length(n.points) != 
        1 || n.points != trunc(n.points) || n.points < 2) 
        stop("'n.points' must be an integer larger than 1")
    switch(x.var, n = {
        if (min.x < 2 || min.x != trunc(min.x)) stop(paste("When x.var=\"n\" the first element of", 
            "range.x.var must be an integer greater than 1"))
        if (max.x != trunc(max.x)) stop(paste("When x.var=\"n\" the second element of", 
            "range.x.var must be an integer greater than", "the first element of range.x.var"))
    }, p = {
        if (min.x < .Machine$double.eps || min.x > 1 - .Machine$double.eps) stop(paste("When x.var=\"p\" the first element of range.x.var", 
            "must be a positive number between 0 and 1"))
        if (max.x > 1 - .Machine$double.eps) stop(paste("When x.var=\"p\" the second element of", 
            "range.x.var must be an positive number between 0 and 1", 
            "and greater than the first element of range.x.var"))
    }, conf.level = {
        if (min.x < .Machine$double.eps || min.x > 1 - .Machine$double.eps) stop(paste("When x.var=\"conf.level\" the first element of range.x.var", 
            "must be a positive number between 0 and 1"))
        if (max.x > 1 - .Machine$double.eps) stop(paste("When x.var=\"conf.level\" the second element of", 
            "range.x.var must be an positive number between 0 and 1", 
            "and greater than the first element of range.x.var"))
    })
    if (length(lcl.rank) != 1 || !is.vector(lcl.rank, mode = "numeric") || 
        !is.finite(lcl.rank) || lcl.rank != trunc(lcl.rank) || 
        lcl.rank < 0) 
        stop("'lcl.rank' must be a non-negative integer")
    if (length(n.plus.one.minus.ucl.rank) != 1 || !is.vector(n.plus.one.minus.ucl.rank, 
        mode = "numeric") || !is.finite(n.plus.one.minus.ucl.rank) || 
        n.plus.one.minus.ucl.rank != trunc(n.plus.one.minus.ucl.rank) || 
        n.plus.one.minus.ucl.rank < 0) 
        stop("'n.plus.one.minus.ucl.rank' must be a non-negative integer")
    if (x.var != "n" && y.var != "n") {
        if (is.null(n) || !is.finite(n) || !is.vector(n, mode = "numeric") || 
            length(n) != 1 || n < 2 || n != trunc(n)) 
            stop("'n' must be an integer greater than 1")
    }
    if (x.var != "p") 
        if (is.null(p) || !is.finite(p) || !is.vector(p, mode = "numeric") || 
            length(p) != 1 || p <= .Machine$double.eps || p >= 
            1 - .Machine$double.eps) {
            stop("'p' must be a scalar between 0 and 1")
        }
    if (x.var != "conf.level" && y.var != "conf.level") 
        if (is.null(conf.level) || !is.finite(conf.level) || 
            !is.vector(conf.level, mode = "numeric") || length(conf.level) != 
            1 || conf.level <= .Machine$double.eps || conf.level >= 
            1 - .Machine$double.eps) {
            stop("'conf.level' must be a scalar between 0 and 1")
        }
    ci.string <- switch(ci.type, two.sided = "(Two-Sided CI)", 
        lower = "(One-Sided Lower CI)", upper = "(One-Sided Upper CI)")
    n.string <- paste("n =", n)
    conf.string <- paste("Confidence Level = ", format(100 * 
        conf.level, digits = digits), "%", sep = "")
    p.string <- paste(format(100 * p, digits = digits), "'th Percentile", 
        sep = "")
    rank.string <- switch(ci.type, two.sided = paste("Rank(LCL) =", 
        lcl.rank, "and [n + 1 - Rank(UCL)] =", n.plus.one.minus.ucl.rank), 
        lower = paste("Rank(LCL) =", lcl.rank), upper = paste("[n + 1 - Rank(UCL)] =", 
            n.plus.one.minus.ucl.rank))
    if (plot.it) 
        gen.gp.list <- checkGraphicsPars(...)$gen.gp.list
    combo <- paste(c(x.var, y.var), collapse = " & ")
    switch(combo, `n & conf.level` = {
        x <- seq(min.x, max.x, by = ceiling((max.x - min.x + 
            1)/n.points))
        if (is.null(xlab)) xlab <- "Sample Size (n)"
        y <- ciNparConfLevel(n = x, p = p, lcl.rank = lcl.rank, 
            n.plus.one.minus.ucl.rank = n.plus.one.minus.ucl.rank, 
            ci.type = ci.type)
        if (is.null(ylab)) ylab <- "Confidence Level"
        line1 <- "Confidence Level vs. Sample Size for Nonparametric"
        line2 <- paste("Confidence Interval for", p.string, ci.string)
        line3 <- rank.string
    }, `conf.level & n` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(xlab)) xlab <- "Confidence Level"
        y <- ciNparN(p = p, lcl.rank = lcl.rank, n.plus.one.minus.ucl.rank = n.plus.one.minus.ucl.rank, 
            ci.type = ci.type, conf.level = x)
        if (is.null(ylab)) ylab <- "Sample Size (n)"
        line1 <- "Sample Size vs. Confidence Level for Nonparametric"
        line2 <- paste("Confidence Interval for", p.string, ci.string)
        line3 <- rank.string
    }, `p & conf.level` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(xlab)) xlab <- "Quantile (p)"
        y <- ciNparConfLevel(n = n, p = x, lcl.rank = lcl.rank, 
            n.plus.one.minus.ucl.rank = n.plus.one.minus.ucl.rank, 
            ci.type = ci.type)
        if (is.null(ylab)) ylab <- "Confidence Level"
        line1 <- "Confidence Level vs. Quantile for Nonparametric"
        line2 <- paste("Confidence Interval with", n.string, 
            ci.string)
        line3 <- rank.string
    }, `p & n` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(xlab)) xlab <- "Quantile (p)"
        y <- ciNparN(p = x, lcl.rank = lcl.rank, n.plus.one.minus.ucl.rank = n.plus.one.minus.ucl.rank, 
            ci.type = ci.type, conf.level = conf.level)
        if (is.null(ylab)) ylab <- "Sample Size (n)"
        line1 <- "Sample Size vs. Quantile for Nonparametric"
        line2 <- paste("Confidence Interval with", conf.string, 
            ci.string)
        line3 <- rank.string
    })
    if (plot.it) {
        if (!add) {
            plot(x, y, type = "n", main = "", sub = "", ..., 
                xlab = xlab, ylab = ylab)
            if (is.null(main)) {
                mtext(text = line1, side = 3, line = 3, cex = 1.25 * 
                  cex.main, font = 2)
                mtext(text = line2, side = 3, line = 1.5, cex = 1.25 * 
                  cex.main, font = 2)
                mtext(text = line3, side = 3, line = 0.25, cex = cex.main)
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
