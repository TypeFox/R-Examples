plotPredIntNormDesign <-
function (x.var = "n", y.var = "half.width", range.x.var = NULL, 
    n = 25, k = 1, n.mean = 1, half.width = 4 * sigma.hat, sigma.hat = 1, 
    method = "Bonferroni", conf.level = 0.95, round.up = FALSE, 
    n.max = 5000, tol = 1e-07, maxiter = 1000, plot.it = TRUE, 
    add = FALSE, n.points = 100, plot.col = "black", plot.lwd = 3 * 
        par("cex"), plot.lty = 1, digits = .Options$digits, cex.main = par("cex"), 
    ..., main = NULL, xlab = NULL, ylab = NULL, type = "l") 
{
    method <- match.arg(method, c("Bonferroni", "exact"))
    x.var <- match.arg(x.var, c("n", "half.width", "k", "sigma.hat", 
        "conf.level"))
    y.var <- match.arg(y.var, c("half.width", "n"))
    if (x.var == y.var) 
        stop("'x.var' and 'y.var' cannot denote the same quantity")
    if (missing(range.x.var)) {
        range.x.var = switch(x.var, n = c(2, 50), half.width = c(2.5 * 
            sigma.hat, 4 * sigma.hat), k = c(1, 20), sigma.hat = c(0.1, 
            2), conf.level = c(0.5, 0.99))
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
    if (x.var != "n" && y.var != "n") {
        if (is.null(n) || !is.finite(n) || !is.vector(n, mode = "numeric") || 
            length(n) != 1 || n < 2 || n != trunc(n)) 
            stop("'n' must be an integer greater than 1")
    }
    if (x.var != "k") {
        if (is.null(k) || !is.finite(k) || !is.vector(k, mode = "numeric") || 
            length(k) != 1 || k < 1 || k != trunc(k)) 
            stop("'k' must be a positive integer")
    }
    switch(x.var, n = {
        if (min.x < 2 || min.x != trunc(min.x)) stop(paste("When x.var=\"n\" the first element of", 
            "range.x.var must be an integer greater than 1"))
        if (max.x != trunc(max.x)) stop(paste("When x.var=\"n\" the second element of", 
            "range.x.var must be an integer greater than", "the first element of range.x.var"))
        if (min.x < k) stop(paste("When x.var=\"n\", 'min.x' must be", 
            "greater than or equal to 'k'"))
    }, half.width = {
        if (min.x < .Machine$double.eps) stop(paste("When x.var=\"half.width\" the first element of", 
            "range.x.var must be a positive number"))
        if (max.x < .Machine$double.eps) stop(paste("When x.var=\"half.width\" the second element of", 
            "range.x.var must be a positive number greater than", 
            "the first element of range.x.var"))
    }, k = {
        if (min.x < 1 || min.x != trunc(min.x)) stop(paste("When x.var=\"k\" the first element of", 
            "range.x.var must be a positive integer"))
        if (max.x != trunc(max.x)) stop(paste("When x.var=\"k\" the second element of", 
            "range.x.var must be an integer greater than", "the first element of range.x.var"))
        if (y.var != "n" && max.x > n) stop(paste("When x.var=\"k\", max.x must be", 
            "less than or equal to n"))
    }, sigma.hat = {
        if (min.x < .Machine$double.eps) stop(paste("When x.var=\"sigma.hat\" the first element of", 
            "range.x.var must be a positive number"))
        if (max.x < .Machine$double.eps) stop(paste("When x.var=\"sigma.hat\" the second element of", 
            "range.x.var must be a positive number greater than", 
            "the first element of range.x.var"))
    }, conf.level = {
        if (min.x < .Machine$double.eps || min.x > 1 - .Machine$double.eps) stop(paste("When x.var=\"conf.level\" the first element of range.x.var", 
            "must be a positive number between 0 and 1"))
        if (max.x > 1 - .Machine$double.eps) stop(paste("When x.var=\"conf.level\" the second element of", 
            "range.x.var must be an positive number between 0 and 1", 
            "and greater than the first element of range.x.var"))
    })
    if (is.null(n.mean) || !is.finite(n.mean) || !is.vector(n.mean, 
        mode = "numeric") || length(n.mean) != 1 || n.mean < 
        1 || n.mean != trunc(n.mean)) 
        stop("'n.mean' must be a positive integer")
    if (x.var != "conf.level") 
        if (is.null(conf.level) || !is.finite(conf.level) || 
            !is.vector(conf.level, mode = "numeric") || length(conf.level) != 
            1 || conf.level <= .Machine$double.eps || conf.level >= 
            1 - .Machine$double.eps) {
            stop("'conf.level' must be a scalar between 0 and 1")
        }
    if (x.var != "sigma.hat" && (is.null(sigma.hat) || !is.finite(sigma.hat) || 
        !is.vector(sigma.hat, mode = "numeric") || length(sigma.hat) != 
        1 || sigma.hat <= .Machine$double.eps)) 
        stop("'sigma.hat' must be a positive scalar")
    if (x.var != "half.width" && y.var != "half.width" && (is.null(half.width) || 
        !is.finite(half.width) || !is.vector(half.width, mode = "numeric") || 
        length(half.width) != 1 || half.width <= .Machine$double.eps)) 
        stop("'half.width' must be a positive scalar")
    if (!is.vector(n.max, mode = "numeric") || length(n.max) != 
        1 || !is.finite(n.max) || n.max != trunc(n.max) || n.max < 
        2) 
        stop("'n.max' must be a positive integer greater than 1")
    if (!is.vector(maxiter, mode = "numeric") || length(maxiter) != 
        1 || !is.finite(maxiter) || maxiter != trunc(maxiter) || 
        maxiter < 2) 
        stop("'maxiter' must be a positive integer greater than 1")
    n.string <- paste("n =", n)
    n.mean.string <- paste("n (Mean) =", n.mean)
    hw.string <- paste("Half-Width =", format(half.width, digits = digits))
    k.string <- paste("k =", k)
    sh.string <- paste("Sigma Hat =", format(sigma.hat, digits = digits))
    conf.string <- paste("Confidence Level =", format(conf.level, 
        digits = digits))
    method.string <- paste("Based on", ifelse(method == "Bonferroni", 
        "Bonferroni", "Exact"), "Method")
    if (n.mean != 1) 
        line3 <- paste(method.string, "and", n.mean.string)
    else line3 <- method.string
    if (plot.it) 
        gen.gp.list <- checkGraphicsPars(...)$gen.gp.list
    combo <- paste(c(x.var, y.var), collapse = " & ")
    switch(combo, `n & half.width` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(xlab)) xlab <- "Sample Size (n)"
        y <- predIntNormHalfWidth(n = x, k = k, n.mean = n.mean, 
            sigma.hat = sigma.hat, method = method, conf.level = conf.level)
        if (is.null(ylab)) ylab <- "Half-Width"
        line1 <- "Half-Width vs. Sample Size for Prediction Interval"
        line2 <- paste(" with ", k.string, ", ", sh.string, ", and ", 
            conf.string, sep = "")
    }, `half.width & n` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(xlab)) xlab <- "Half-Width"
        y <- predIntNormN(half.width = x, k = k, n.mean = n.mean, 
            sigma.hat = sigma.hat, method = method, conf.level = conf.level, 
            round.up = round.up, n.max = n.max, tol = tol, maxiter = maxiter)
        if (is.null(ylab)) ylab <- "Sample Size (n)"
        line1 <- "Sample Size vs. Half-Width for Prediction Interval"
        line2 <- paste(" with ", k.string, ", ", sh.string, ", and ", 
            conf.string, sep = "")
    }, `k & half.width` = {
        x <- seq(min.x, max.x, by = ceiling((max.x - min.x + 
            1)/n.points))
        if (is.null(xlab)) xlab <- "# Future Observations (k)"
        y <- predIntNormHalfWidth(n = n, k = x, n.mean = n.mean, 
            sigma.hat = sigma.hat, method = method, conf.level = conf.level)
        if (is.null(ylab)) ylab <- "Half-Width"
        line1 <- "Half-Width vs. # Future Observations for Prediction Interval"
        line2 <- paste(" with ", n.string, ", ", sh.string, ", and ", 
            conf.string, sep = "")
    }, `k & n` = {
        x <- seq(min.x, max.x, by = ceiling((max.x - min.x + 
            1)/n.points))
        if (is.null(xlab)) xlab <- "# Future Observations (k)"
        y <- predIntNormN(half.width = half.width, k = x, n.mean = n.mean, 
            sigma.hat = sigma.hat, method = method, conf.level = conf.level, 
            round.up = round.up, n.max = n.max, tol = tol, maxiter = maxiter)
        if (is.null(ylab)) ylab <- "Sample Size (n)"
        line1 <- "Sample Size vs. # Future Observations for Prediction Interval"
        line2 <- paste(" with ", hw.string, ", ", sh.string, 
            ", and ", conf.string, sep = "")
    }, `sigma.hat & half.width` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(xlab)) xlab <- "Sigma Hat"
        y <- predIntNormHalfWidth(n = n, k = k, n.mean = n.mean, 
            sigma.hat = x, method = method, conf.level = conf.level)
        if (is.null(ylab)) ylab <- "Half-Width"
        line1 <- "Half-Width vs. Sigma Hat for Prediction Interval"
        line2 <- paste(" with ", n.string, ", ", k.string, ", and ", 
            conf.string, sep = "")
    }, `sigma.hat & n` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(xlab)) xlab <- "Sigma Hat"
        y <- predIntNormN(half.width = half.width, k = k, n.mean = n.mean, 
            sigma.hat = x, method = method, conf.level = conf.level, 
            round.up = round.up, n.max = n.max, tol = tol, maxiter = maxiter)
        if (is.null(ylab)) ylab <- "Sample Size (n)"
        line1 <- "Sample Size vs. Sigma Hat for Prediction Interval"
        line2 <- paste(" with ", k.string, ", ", hw.string, ", and ", 
            conf.string, sep = "")
    }, `conf.level & half.width` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(xlab)) xlab <- "Confidence Level"
        y <- predIntNormHalfWidth(n = n, k = k, n.mean = n.mean, 
            sigma.hat = sigma.hat, method = method, conf.level = x)
        if (is.null(ylab)) ylab <- "Half-Width"
        line1 <- "Half-Width vs. Confidence Level for Prediction Interval"
        line2 <- paste(" with ", n.string, ", ", k.string, ", and ", 
            sh.string, sep = "")
    }, `conf.level & n` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(xlab)) xlab <- "Confidence Level"
        y <- predIntNormN(half.width = half.width, k = k, n.mean = n.mean, 
            sigma.hat = sigma.hat, method = method, conf.level = x, 
            round.up = round.up, n.max = n.max, tol = tol, maxiter = maxiter)
        if (is.null(ylab)) ylab <- "Sample Size (n)"
        line1 <- "Sample Size vs. Confidence Level for Prediction Interval"
        line2 <- paste(" with ", k.string, ", ", hw.string, ", and ", 
            sh.string, sep = "")
    })
    if (plot.it) {
        if (!add) {
            plot(x, y, type = "n", main = "", sub = "", ..., 
                xlab = xlab, ylab = ylab)
            if (is.null(main)) {
                mtext(text = line1, side = 3, line = 3, cex = 1.25 * 
                  cex.main)
                mtext(text = line2, side = 3, line = 1.75, cex = 1.25 * 
                  cex.main)
                mtext(text = line3, side = 3, line = 0.5, cex = cex.main)
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
