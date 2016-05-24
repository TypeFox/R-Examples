plotTolIntNormDesign <-
function (x.var = "n", y.var = "half.width", range.x.var = NULL, 
    n = 25, half.width = ifelse(x.var == "sigma.hat", 3 * max.x, 
        3 * sigma.hat), sigma.hat = 1, coverage = 0.95, conf.level = 0.95, 
    cov.type = "content", round.up = FALSE, n.max = 5000, tol = 1e-07, 
    maxiter = 1000, plot.it = TRUE, add = FALSE, n.points = 100, 
    plot.col = 1, plot.lwd = 3 * par("cex"), plot.lty = 1, digits = .Options$digits, 
    ..., main = NULL, xlab = NULL, ylab = NULL, type = "l") 
{
    x.var <- match.arg(x.var, c("n", "half.width", "sigma.hat", 
        "coverage", "conf.level"))
    y.var <- match.arg(y.var, c("half.width", "n"))
    cov.type <- match.arg(cov.type, c("content", "expectation"))
    if (x.var == y.var) 
        stop("'x.var' and 'y.var' cannot denote the same quantity")
    if (missing(range.x.var)) {
        range.x.var <- switch(x.var, n = c(2, 50), half.width = c(2.5 * 
            sigma.hat, 4 * sigma.hat), sigma.hat = c(0.1, 2), 
            coverage = c(0.5, 0.99), conf.level = c(0.5, 0.99))
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
    if (x.var != "sigma.hat") {
        if (length(sigma.hat) != 1 || !is.vector(sigma.hat, mode = "numeric") || 
            !is.finite(sigma.hat) || sigma.hat < .Machine$double.eps) 
            stop("'sigma.hat' must be a positive numeric scalar")
    }
    switch(x.var, n = {
        if (min.x < 2 || min.x != trunc(min.x)) stop(paste("When x.var=\"n\" the first element of", 
            "range.x.var must be an integer greater than 1"))
        if (max.x != trunc(max.x)) stop(paste("When x.var=\"n\" the second element of", 
            "range.x.var must be an integer greater than", "the first element of range.x.var"))
    }, half.width = {
        if (min.x < .Machine$double.eps) stop(paste("When x.var=\"half.width\" the first element of", 
            "range.x.var must be a positive number"))
        if (max.x < .Machine$double.eps) stop(paste("When x.var=\"half.width\" the second element of", 
            "range.x.var must be a positive number greater than", 
            "the first element of range.x.var"))
    }, sigma.hat = {
        if (min.x < .Machine$double.eps) stop(paste("When x.var=\"sigma.hat\" the first element of", 
            "range.x.var must be a positive number"))
        if (max.x < .Machine$double.eps) stop(paste("When x.var=\"sigma.hat\" the second element of", 
            "range.x.var must be a positive number greater than", 
            "the first element of range.x.var"))
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
    if (x.var != "n" && y.var != "n") {
        if (is.null(n) || !is.finite(n) || !is.vector(n, mode = "numeric") || 
            length(n) != 1 || n < 2 || n != trunc(n)) 
            stop("'n' must be an integer greater than 1")
    }
    if (x.var != "half.width" && y.var != "half.width") {
        if (length(half.width) != 1 || !is.finite(half.width) || 
            !is.vector(half.width, mode = "numeric") || half.width <= 
            .Machine$double.eps) 
            stop("'half.width' must be a positive scalar")
    }
    if (x.var != "coverage") 
        if (is.null(coverage) || !is.finite(coverage) || !is.vector(coverage, 
            mode = "numeric") || length(coverage) != 1 || coverage <= 
            .Machine$double.eps || coverage >= 1 - .Machine$double.eps) {
            stop("'coverage' must be a scalar between 0 and 1")
        }
    if (x.var != "conf.level") 
        if (is.null(conf.level) || !is.finite(conf.level) || 
            !is.vector(conf.level, mode = "numeric") || length(conf.level) != 
            1 || conf.level <= .Machine$double.eps || conf.level >= 
            1 - .Machine$double.eps) {
            stop("'conf.level' must be a scalar between 0 and 1")
        }
    if (!is.vector(n.max, mode = "numeric") || length(n.max) != 
        1 || !is.finite(n.max) || n.max != trunc(n.max) || n.max < 
        2) 
        stop("'n.max' must be a positive integer greater than 1")
    if (!is.vector(maxiter, mode = "numeric") || length(maxiter) != 
        1 || !is.finite(maxiter) || maxiter != trunc(maxiter) || 
        maxiter < 2) 
        stop("'maxiter' must be a positive integer greater than 1")
    cov.type.string <- ifelse(cov.type == "content", "B-Content Tolerance Interval", 
        "B-Expectation Tolerance Interval")
    n.string <- paste("n =", n)
    hw.string <- paste("Half-Width =", format(half.width, digits = digits))
    sh.string <- paste("Sigma Hat =", format(sigma.hat, digits = digits))
    cov.string <- paste("Coverage = ", format(100 * coverage, 
        digits = digits), "%", sep = "")
    conf.string <- paste("Confidence Level = ", format(100 * 
        conf.level, digits = digits), "%", sep = "")
    if (plot.it) 
        gen.gp.list <- checkGraphicsPars(...)$gen.gp.list
    combo <- paste(c(x.var, y.var), collapse = " & ")
    switch(combo, `n & half.width` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(xlab)) xlab <- "Sample Size (n)"
        y <- tolIntNormHalfWidth(n = x, sigma.hat = sigma.hat, 
            coverage = coverage, cov.type = cov.type, conf.level = conf.level)
        if (is.null(ylab)) ylab <- "Half-Width"
        line1 <- paste("Half-Width vs. Sample Size for", cov.type.string)
        if (cov.type == "content") line2 <- paste(" with ", sh.string, 
            ", ", cov.string, ", and ", conf.string, sep = "") else line2 <- paste(" with ", 
            sh.string, " and ", cov.string, sep = "")
    }, `half.width & n` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(xlab)) xlab <- "Half-Width"
        y <- tolIntNormN(half.width = x, sigma.hat = sigma.hat, 
            coverage = coverage, cov.type = cov.type, conf.level = conf.level, 
            round.up = round.up, n.max = n.max, tol = tol, maxiter = maxiter)
        if (is.null(ylab)) ylab <- "Sample Size (n)"
        line1 <- paste("Sample Size vs. Half-Width for", cov.type.string)
        if (cov.type == "content") line2 <- paste(" with ", sh.string, 
            ", ", cov.string, ", and ", conf.string, sep = "") else line2 <- paste(" with ", 
            sh.string, " and ", cov.string, sep = "")
    }, `sigma.hat & half.width` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(xlab)) xlab <- "Sigma Hat"
        y <- tolIntNormHalfWidth(n = n, sigma.hat = x, coverage = coverage, 
            cov.type = cov.type, conf.level = conf.level)
        if (is.null(ylab)) ylab <- "Half-Width"
        line1 <- paste("Half-Width vs. Sigma Hat for", cov.type.string)
        if (cov.type == "content") line2 <- paste(" with ", n.string, 
            ", ", cov.string, ", and ", conf.string, sep = "") else line2 <- paste(" with ", 
            n.string, " and ", cov.string, sep = "")
    }, `sigma.hat & n` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(xlab)) xlab <- "Sigma Hat"
        y <- tolIntNormN(half.width = half.width, sigma.hat = x, 
            coverage = coverage, cov.type = cov.type, conf.level = conf.level, 
            round.up = round.up, n.max = n.max, tol = tol, maxiter = maxiter)
        if (is.null(ylab)) ylab <- "Sample Size (n)"
        line1 <- paste("Sample Size vs. Sigma Hat for", cov.type.string)
        if (cov.type == "content") line2 <- paste(" with ", hw.string, 
            ", ", cov.string, ", and ", conf.string, sep = "") else line2 <- paste(" with ", 
            hw.string, " and ", cov.string, sep = "")
    }, `coverage & half.width` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(xlab)) xlab <- "Coverage"
        y <- tolIntNormHalfWidth(n = n, sigma.hat = sigma.hat, 
            coverage = x, cov.type = cov.type, conf.level = conf.level)
        if (is.null(ylab)) ylab <- "Half-Width"
        line1 <- paste("Half-Width vs. Coverage for", cov.type.string)
        if (cov.type == "content") line2 <- paste("with ", n.string, 
            ", ", sh.string, ", and ", conf.string, sep = "") else line2 <- paste("with ", 
            n.string, " and ", sh.string, sep = "")
    }, `coverage & n` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(xlab)) xlab <- "Coverage"
        y <- tolIntNormN(half.width = half.width, sigma.hat = sigma.hat, 
            coverage = x, cov.type = cov.type, conf.level = conf.level, 
            round.up = round.up, n.max = n.max, tol = tol, maxiter = maxiter)
        if (is.null(ylab)) ylab <- "Sample Size (n)"
        line1 <- paste("Sample Size vs. Coverage for", cov.type.string)
        if (cov.type == "content") line2 <- paste("with ", hw.string, 
            ", ", sh.string, ", and ", conf.string, sep = "") else line2 <- paste("with ", 
            hw.string, " and ", sh.string, sep = "")
    }, `conf.level & half.width` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(xlab)) xlab <- "Confidence Level"
        y <- tolIntNormHalfWidth(n = n, sigma.hat = sigma.hat, 
            coverage = coverage, cov.type = cov.type, conf.level = x)
        if (is.null(ylab)) ylab <- "Half-Width"
        line1 <- paste("Half-Width vs. Confidence Level for", 
            cov.type.string)
        line2 <- paste("with ", n.string, ", ", sh.string, ", and ", 
            cov.string, sep = "")
    }, `conf.level & n` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(xlab)) xlab <- "Confidence Level"
        y <- tolIntNormN(half.width = half.width, sigma.hat = sigma.hat, 
            coverage = coverage, cov.type = cov.type, conf.level = x, 
            round.up = round.up, n.max = n.max, tol = tol, maxiter = maxiter)
        if (is.null(ylab)) ylab <- "Sample Size (n)"
        line1 <- paste("Sample Size vs. Confidence Level for", 
            cov.type.string)
        line2 <- paste("with ", hw.string, ", ", sh.string, ", and ", 
            cov.string, sep = "")
    })
    if (plot.it) {
        if (!add) {
            plot(x, y, type = "n", main = "", sub = "", ..., 
                xlab = xlab, ylab = ylab)
            if (is.null(main)) 
                main <- paste(line1, line2, sep = "\n")
            arg.list <- c(list(main = main), gen.gp.list)
            do.call("title", arg.list)
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
