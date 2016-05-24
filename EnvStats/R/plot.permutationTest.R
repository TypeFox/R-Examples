plot.permutationTest <-
function (x, hist.col = "cyan", stat.col = "black", stat.lwd = 3 * 
    par("cex"), stat.lty = 1, cex.main = par("cex"), digits = .Options$digits, 
    main = NULL, xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, 
    ...) 
{
    exact <- x$exact
    subset.expression <- x$subset.expression
    is.null.subset.expression <- is.null(subset.expression)
    if (is.null(main)) {
        if (exact) {
            num.top.lines <- ifelse(is.null.subset.expression, 
                4, 5)
        }
        else {
            num.top.lines <- ifelse(is.null.subset.expression, 
                5, 6)
        }
        o.mar <- par(mar = c(5, 4, num.top.lines, 4) + 0.1)
        on.exit(par(o.mar))
    }
    if (is.null(xlab)) 
        xlab <- names(x$statistic)
    if (is.null(ylab)) 
        ylab <- "Relative Frequency"
    stat <- x$statistic
    stat.dist <- x$stat.dist
    probs.stat.dist <- x$probs.stat.dist
    if (is.null(probs.stat.dist)) {
        hist.list <- do.call("hist", c(list(x = stat.dist, plot = FALSE), 
            list(...)))
        if (is.null(xlim)) {
            xlim <- range(stat.dist, stat, hist.list$breaks)
        }
        do.call("hist", c(list(x = stat.dist, prob = TRUE, col = hist.col), 
            list(...), list(xlab = xlab, ylab = ylab, xlim = xlim, 
                main = "")))
    }
    else {
        nx <- length(stat.dist)
        min.diff <- min(diff(stat.dist))
        con.upper.bound <- (min.diff/2)
        con <- con.upper.bound * (0.9 + ((nx - 2)/nx) * 0.1)
        xleft <- stat.dist - con
        xright <- stat.dist + con
        ybottom <- rep(0, nx)
        if (is.null(xlim)) {
            xlim <- range(stat.dist, stat, xleft, xright)
        }
        if (is.null(ylim)) {
            ylim <- c(0, max(probs.stat.dist))
        }
        plot(xleft, probs.stat.dist, type = "n", ..., xlim = xlim, 
            ylim = ylim, xlab = xlab, ylab = ylab)
        rect(xleft = xleft, ybottom = ybottom, xright = xright, 
            ytop = probs.stat.dist, col = hist.col, border = stat.col, 
            ...)
    }
    abline(v = x$statistic, col = stat.col, lwd = stat.lwd, lty = stat.lty)
    if (is.null(main)) {
        nv <- x$null.value
        nnv <- names(nv)
        method.vec <- string.break.line(x$method)[[1]]
        data.name <- paste(x$data.name, collapse = " and ")
        parent.of.data <- x$parent.of.data
        if (!is.null(parent.of.data)) 
            data.name <- paste(data.name, "in", parent.of.data)
        main1 <- paste("Results of", method.vec[1], "for")
        main2 <- data.name
        if (!is.null.subset.expression) {
            main3 <- paste("(Subset with ", subset.expression, 
                ")", sep = "")
        }
        alternative <- x$alternative
        string <- switch(alternative, two.sided = "Two-Sided", 
            less = "One-Sided Lower", greater = "One-Sided Upper")
        string <- paste(string, "P-Value =", signif(x$p.value, 
            digits = digits))
        main.H0 <- paste("H0:", paste(format(nnv, justify = "left"), 
            format(nv, digits = digits), sep = " = "))
        main.H0 <- paste(main.H0, "; ", string)
        if (exact) {
            if (is.null.subset.expression) {
                mtext(main1, side = 3, line = 3, cex = cex.main)
                mtext(main2, side = 3, line = 2, cex = cex.main)
                mtext(main.H0, side = 3, line = 1, cex = cex.main)
            }
            else {
                mtext(main1, side = 3, line = 4, cex = cex.main)
                mtext(main2, side = 3, line = 3, cex = cex.main)
                mtext(main3, side = 3, line = 2, cex = cex.main)
                mtext(main.H0, side = 3, line = 1, cex = cex.main)
            }
        }
        else {
            method.vec <- method.vec[-1]
            nchar.vec <- nchar(method.vec)
            main4 <- paste(substring(method.vec, 34, nchar.vec), 
                collapse = " ")
            if (is.null.subset.expression) {
                mtext(main1, side = 3, line = 4, cex = cex.main)
                mtext(main2, side = 3, line = 3, cex = cex.main)
                mtext(main4, side = 3, line = 2, cex = cex.main)
                mtext(main.H0, side = 3, line = 1, cex = cex.main)
            }
            else {
                mtext(main1, side = 3, line = 5, cex = cex.main)
                mtext(main2, side = 3, line = 4, cex = cex.main)
                mtext(main3, side = 3, line = 3, cex = cex.main)
                mtext(main4, side = 3, line = 2, cex = cex.main)
                mtext(main.H0, side = 3, line = 1, cex = cex.main)
            }
        }
    }
    else title(main = main)
    invisible(x)
}
