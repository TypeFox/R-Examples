plotTolIntNparDesign <-
function (x.var = "n", y.var = "conf.level", range.x.var = NULL, 
    n = 25, coverage = 0.95, conf.level = 0.95, ti.type = "two.sided", 
    cov.type = "content", ltl.rank = ifelse(ti.type == "upper", 
        0, 1), n.plus.one.minus.utl.rank = ifelse(ti.type == 
        "lower", 0, 1), plot.it = TRUE, add = FALSE, n.points = 100, 
    plot.col = "black", plot.lwd = 3 * par("cex"), plot.lty = 1, 
    digits = .Options$digits, cex.main = par("cex"), ..., main = NULL, 
    xlab = NULL, ylab = NULL, type = "l") 
{
    x.var <- match.arg(x.var, c("n", "coverage", "conf.level"))
    y.var <- match.arg(y.var, c("conf.level", "n", "coverage"))
    ti.type <- match.arg(ti.type, c("two.sided", "lower", "upper"))
    cov.type <- match.arg(cov.type, c("content", "expectation"))
    if (x.var == y.var) 
        stop("'x.var' and 'y.var' cannot denote the same quantity")
    if (missing(range.x.var)) {
        range.x.var <- switch(x.var, n = c(2, 50), coverage = c(0.5, 
            0.99), conf.level = c(0.5, 0.99))
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
    }, coverage = {
        if (min.x < .Machine$double.eps || min.x > 1 - .Machine$double.eps) stop(paste("When x.var=\"coverage\" the first element of range.x.var", 
            "must be a positive number between 0 and 1"))
        if (max.x > 1 - .Machine$double.eps) stop(paste("When x.var=\"coverage\" the second element of", 
            "range.x.var must be an positive number between 0 and 1", 
            "and greater than the first element of range.x.var"))
    }, conf.level = {
        if (min.x < .Machine$double.eps || min.x > 1 - .Machine$double.eps) stop(paste("When x.var=\"conf.level\" the first element of range.x.var", 
            "must be a positive number between 0 and 1"))
        if (max.x > 1 - .Machine$double.eps) stop(paste("When x.var=\"conf.level\" the second element of", 
            "range.x.var must be an positive number between 0 and 1", 
            "and greater than the first element of range.x.var"))
        if (cov.type != "content") {
            cov.type <- "content"
            warning(paste("When x.var=\"conf.level\" cov.type", 
                "is automatically set to \"content\""))
        }
    })
    if (length(ltl.rank) != 1 || !is.vector(ltl.rank, mode = "numeric") || 
        !is.finite(ltl.rank) || ltl.rank != trunc(ltl.rank) || 
        ltl.rank < 0) 
        stop("'ltl.rank' must be a non-negative integer")
    if (length(n.plus.one.minus.utl.rank) != 1 || !is.vector(n.plus.one.minus.utl.rank, 
        mode = "numeric") || !is.finite(n.plus.one.minus.utl.rank) || 
        n.plus.one.minus.utl.rank != trunc(n.plus.one.minus.utl.rank) || 
        n.plus.one.minus.utl.rank < 0) 
        stop("'n.plus.one.minus.utl.rank' must be a non-negative integer")
    if (y.var == "conf.level" && cov.type != "content") {
        cov.type <- "content"
        warning(paste("When y.var=\"conf.level\" cov.type", "is automatically set to \"content\""))
    }
    if (x.var != "n" && y.var != "n") {
        if (is.null(n) || !is.finite(n) || !is.vector(n, mode = "numeric") || 
            length(n) != 1 || n < 2 || n != trunc(n)) 
            stop("'n' must be an integer greater than 1")
    }
    if (x.var != "coverage" && y.var != "coverage") 
        if (is.null(coverage) || !is.finite(coverage) || !is.vector(coverage, 
            mode = "numeric") || length(coverage) != 1 || coverage <= 
            .Machine$double.eps || coverage >= 1 - .Machine$double.eps) {
            stop("'coverage' must be a scalar between 0 and 1")
        }
    if (x.var != "conf.level" && y.var != "conf.level") 
        if (is.null(conf.level) || !is.finite(conf.level) || 
            !is.vector(conf.level, mode = "numeric") || length(conf.level) != 
            1 || conf.level <= .Machine$double.eps || conf.level >= 
            1 - .Machine$double.eps) {
            stop("'conf.level' must be a scalar between 0 and 1")
        }
    ti.string <- switch(ti.type, two.sided = "(Two-Sided TI)", 
        lower = "(One-Sided Lower TI)", upper = "(One-Sided Upper TI)")
    cov.type.string <- ifelse(cov.type == "content", "Nonparametric B-Content Tolerance Interval", 
        "Nonparametric B-Expectation Tolerance Interval")
    n.string <- paste("n =", n)
    conf.string <- paste("Confidence Level = ", format(100 * 
        conf.level, digits = digits), "%", sep = "")
    cov.string <- paste("Coverage = ", format(100 * coverage, 
        digits = digits), "%", sep = "")
    rank.string <- switch(ti.type, two.sided = paste("Rank(LTL) =", 
        ltl.rank, "and [n + 1 - Rank(UTL)] =", n.plus.one.minus.utl.rank), 
        lower = paste("Rank(LTL) =", ltl.rank), upper = paste("[n + 1 - Rank(UTL)] =", 
            n.plus.one.minus.utl.rank))
    if (plot.it) 
        gen.gp.list <- checkGraphicsPars(...)$gen.gp.list
    combo <- paste(c(x.var, y.var), collapse = " & ")
    switch(combo, `n & conf.level` = {
        x <- seq(min.x, max.x, by = ceiling((max.x - min.x + 
            1)/n.points))
        if (is.null(xlab)) xlab <- "Sample Size (n)"
        y <- tolIntNparConfLevel(n = x, coverage = coverage, 
            ltl.rank = ltl.rank, n.plus.one.minus.utl.rank = n.plus.one.minus.utl.rank, 
            ti.type = ti.type)
        if (is.null(ylab)) ylab <- "Confidence Level"
        line1 <- "Confidence Level vs. Sample Size for"
        line2 <- paste(cov.type.string, "with", cov.string, ti.string)
        line3 <- rank.string
    }, `n & coverage` = {
        x <- seq(min.x, max.x, by = ceiling((max.x - min.x + 
            1)/n.points))
        if (is.null(xlab)) xlab <- "Sample Size (n)"
        y <- tolIntNparCoverage(n = x, conf.level = conf.level, 
            cov.type = cov.type, ltl.rank = ltl.rank, n.plus.one.minus.utl.rank = n.plus.one.minus.utl.rank, 
            ti.type = ti.type)
        if (is.null(ylab)) ylab <- "Coverage"
        line1 <- "Coverage vs. Sample Size for"
        line2 <- ifelse(cov.type == "content", paste(cov.type.string, 
            "with", conf.string, ti.string), paste(cov.type.string, 
            ti.string))
        line3 <- rank.string
    }, `conf.level & n` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(xlab)) xlab <- "Confidence Level"
        y <- tolIntNparN(ltl.rank = ltl.rank, n.plus.one.minus.utl.rank = n.plus.one.minus.utl.rank, 
            coverage = coverage, cov.type = cov.type, ti.type = ti.type, 
            conf.level = x)
        if (is.null(ylab)) ylab <- "Sample Size (n)"
        line1 <- "Sample Size vs. Confidence Level for"
        line2 <- paste(cov.type.string, "with", cov.string, ti.string)
        line3 <- rank.string
    }, `conf.level & coverage` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(xlab)) xlab <- "Confidence Level"
        y <- tolIntNparCoverage(n = n, conf.level = x, cov.type = cov.type, 
            ltl.rank = ltl.rank, n.plus.one.minus.utl.rank = n.plus.one.minus.utl.rank, 
            ti.type = ti.type)
        if (is.null(ylab)) ylab <- "Coverage"
        line1 <- "Coverage vs. Confidence Level for"
        line2 <- paste(cov.type.string, "with", n.string, ti.string)
        line3 <- rank.string
    }, `coverage & conf.level` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(xlab)) xlab <- "Coverage"
        y <- tolIntNparConfLevel(n = n, coverage = x, ltl.rank = ltl.rank, 
            n.plus.one.minus.utl.rank = n.plus.one.minus.utl.rank, 
            ti.type = ti.type)
        if (is.null(ylab)) ylab <- "Confidence Level"
        line1 <- "Confidence Level vs. Coverage for"
        line2 <- paste(cov.type.string, "with", n.string, ti.string)
        line3 <- rank.string
    }, `coverage & n` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(xlab)) xlab <- "Coverage"
        y <- tolIntNparN(ltl.rank = ltl.rank, n.plus.one.minus.utl.rank = n.plus.one.minus.utl.rank, 
            coverage = x, cov.type = cov.type, ti.type = ti.type, 
            conf.level = conf.level)
        if (is.null(ylab)) ylab <- "Sample Size (n)"
        line1 <- "Sample Size vs. Coverage for"
        line2 <- ifelse(cov.type == "content", paste(cov.type.string, 
            "with", conf.string, ti.string), paste(cov.type.string, 
            ti.string))
        line3 <- rank.string
    })
    if (plot.it) {
        if (!add) {
            plot(x, y, type = "n", main = "", sub = "", ..., 
                xlab = xlab, ylab = ylab)
            if (is.null(main)) {
                mtext(text = line1, side = 3, line = 3, cex = 1.2 * 
                  cex.main)
                mtext(text = line2, side = 3, line = 1.75, cex = 1.2 * 
                  cex.main)
                mtext(text = line3, side = 3, line = 0.5, cex = 0.9 * 
                  cex.main)
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
