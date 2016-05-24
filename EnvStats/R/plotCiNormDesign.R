plotCiNormDesign <-
function (x.var = "n", y.var = "half.width", range.x.var = NULL, 
    n.or.n1 = 25, n2 = n.or.n1, half.width = sigma.hat/2, sigma.hat = 1, 
    conf.level = 0.95, sample.type = ifelse(missing(n2), "one.sample", 
        "two.sample"), round.up = FALSE, n.max = 5000, tol = 1e-07, 
    maxiter = 1000, plot.it = TRUE, add = FALSE, n.points = 100, 
    plot.col = "black", plot.lwd = 3 * par("cex"), plot.lty = 1, 
    digits = .Options$digits, main = NULL, xlab = NULL, ylab = NULL, 
    type = "l", ...) 
{
    x.var <- match.arg(x.var, c("n", "half.width", "sigma.hat", 
        "conf.level"))
    y.var <- match.arg(y.var, c("half.width", "n"))
    if (x.var == y.var) 
        stop("'x.var' and 'y.var' cannot denote the same quantity")
    if (missing(range.x.var)) {
        range.x.var <- switch(x.var, n = c(2, 50), half.width = c(0.1/sigma.hat, 
            2/sigma.hat), sigma.hat = c(0.1, 2), conf.level = c(0.5, 
            0.99))
    }
    if (is.null(range.x.var) || !all(is.finite(range.x.var)) || 
        !is.vector(range.x.var, mode = "numeric") || length(range.x.var) != 
        2) 
        stop(paste("'range.x.var' must be a numeric vector of length 2", 
            "with no missing (NA), infinite(-Inf, Inf), or undefined(NaN) values"))
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
    }, conf.level = {
        if (min.x < .Machine$double.eps || min.x > 1 - .Machine$double.eps) stop(paste("When x.var=\"conf.level\" the first element of range.x.var", 
            "must be a positive number between 0 and 1"))
        if (max.x > 1 - .Machine$double.eps) stop(paste("When x.var=\"conf.level\" the second element of", 
            "range.x.var must be an positive number between 0 and 1", 
            "and greater than the first element of range.x.var"))
    })
    if (x.var != "conf.level") 
        if (is.null(conf.level) || !is.finite(conf.level) || 
            !is.vector(conf.level, mode = "numeric") || length(conf.level) != 
            1 || conf.level <= .Machine$double.eps || conf.level >= 
            1 - .Machine$double.eps) {
            stop("'conf.level' must be a scalar between 0 and 1")
        }
    sample.type <- match.arg(sample.type, c("one.sample", "two.sample"))
    n2.constrained <- !missing(n2) && !is.null(n2)
    if (x.var != "n" && y.var != "n") {
        if (is.null(n.or.n1) || !is.finite(n.or.n1) || !is.vector(n.or.n1, 
            mode = "numeric") || length(n.or.n1) != 1 || n.or.n1 < 
            2 || n.or.n1 != trunc(n.or.n1)) 
            stop("'n.or.n1' must be an integer greater than 1")
        if (sample.type == "two.sample" && (is.null(n2) || !is.finite(n2) || 
            !is.vector(n2, mode = "numeric") || length(n2) != 
            1 || n2 < 2 || n2 != trunc(n2))) 
            stop("'n2' must be an integer greater than 1")
    }
    else if (sample.type == "two.sample" && n2.constrained && 
        (is.null(n2) || !is.finite(n2) || !is.vector(n2, mode = "numeric") || 
            length(n2) != 1 || n2 < 2 || n2 != trunc(n2))) 
        stop("'n2' must be an integer greater than 1")
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
    type.string <- ifelse(sample.type == "one.sample", "Mean,", 
        "(Mean 1 - Mean 2),")
    if (plot.it) 
        gen.gp.list <- checkGraphicsPars(...)$gen.gp.list
    combo <- paste(sort(c(x.var, y.var)), collapse = " & ")
    switch(combo, `conf.level & half.width` = {
        x <- seq(min.x, max.x, length = n.points)
        y <- ciNormHalfWidth(n.or.n1 = n.or.n1, n2 = n2, sigma.hat = sigma.hat, 
            conf.level = x, sample.type = sample.type)
        if (is.null(ylab)) ylab <- "Half-Width"
        if (is.null(xlab)) xlab <- "Confidence Level"
        if (sample.type == "two.sample") {
            n.string <- paste("n1 = ", n.or.n1, ", n2 = ", n2, 
                ",", sep = "")
        } else {
            n.string <- paste("n =", n.or.n1)
        }
        if (is.null(main)) main <- paste("Half-Width vs. Confidence Level for ", 
            "Confidence Interval for\n", type.string, " with ", 
            n.string, " and Estimated SD = ", format(sigma.hat, 
                digits = digits), sep = "")
    }, `conf.level & n` = {
        x <- seq(min.x, max.x, length = n.points)
        if (sample.type == "two.sample") {
            if (n2.constrained) {
                if (is.null(ylab)) ylab <- paste("Sample Size (n1), With n2 =", 
                  n2)
                y <- ciNormN(half.width = half.width, sigma.hat = sigma.hat, 
                  conf.level = x, sample.type = "two.sample", 
                  n2 = n2, round.up = round.up, n.max = n.max, 
                  tol = tol, maxiter = maxiter)$n1
            } else {
                if (is.null(ylab)) ylab <- "Sample Size (n1 and n2)"
                y <- ciNormN(half.width = half.width, sigma.hat = sigma.hat, 
                  conf.level = x, sample.type = "two.sample", 
                  round.up = round.up, n.max = n.max, tol = tol, 
                  maxiter = maxiter)
            }
        } else {
            if (is.null(ylab)) ylab <- "Sample Size (n)"
            y <- ciNormN(half.width = half.width, sigma.hat = sigma.hat, 
                conf.level = x, sample.type = "one.sample", round.up = round.up, 
                n.max = n.max, tol = tol, maxiter = maxiter)
        }
        if (is.null(xlab)) xlab <- "Confidence Level"
        if (is.null(main)) main <- paste("Sample Size vs. Confidence Level for ", 
            "Confidence Interval for\n", type.string, " with Estimated SD = ", 
            format(sigma.hat, digits = digits), " and Half-Width = ", 
            format(half.width, digits = digits), sep = "")
    }, `n & sigma.hat` = {
        x <- seq(min.x, max.x, length = n.points)
        if (sample.type == "two.sample") {
            if (n2.constrained) {
                if (is.null(ylab)) ylab <- paste("Sample Size (n1), With n2 =", 
                  n2)
                y <- ciNormN(half.width = half.width, sigma.hat = x, 
                  conf.level = conf.level, sample.type = "two.sample", 
                  n2 = n2, round.up = round.up, n.max = n.max, 
                  tol = tol, maxiter = maxiter)$n1
            } else {
                if (is.null(ylab)) ylab <- "Sample Size (n1 and n2)"
                y <- ciNormN(half.width = half.width, sigma.hat = x, 
                  conf.level = conf.level, sample.type = "two.sample", 
                  round.up = round.up, n.max = n.max, tol = tol, 
                  maxiter = maxiter)
            }
        } else {
            if (is.null(ylab)) ylab <- "Sample Size (n)"
            y <- ciNormN(half.width = half.width, sigma.hat = x, 
                conf.level = conf.level, sample.type = "one.sample", 
                round.up = round.up, n.max = n.max, tol = tol, 
                maxiter = maxiter)
        }
        if (is.null(xlab)) xlab <- "Estimated SD"
        if (is.null(main)) main <- paste("Sample Size vs. Estimated SD for ", 
            "Confidence Interval for\n", type.string, " with Half-Width = ", 
            format(half.width, digits = digits), " and Confidence Level = ", 
            format(conf.level, digits = digits), sep = "")
    }, `half.width & sigma.hat` = {
        x <- seq(min.x, max.x, length = n.points)
        y <- ciNormHalfWidth(n.or.n1 = n.or.n1, n2 = n2, sigma.hat = x, 
            conf.level = conf.level, sample.type = sample.type)
        if (is.null(xlab)) xlab <- "Estimated SD"
        if (is.null(ylab)) ylab <- "Half-Width"
        if (sample.type == "two.sample") {
            n.string <- paste("n1 = ", n.or.n1, ", n2 = ", n2, 
                ",", sep = "")
        } else {
            n.string <- paste("n =", n.or.n1)
        }
        if (is.null(main)) main <- paste("Half-Width vs. Estimated SD for ", 
            "Confidence Interval for\n", type.string, " with ", 
            n.string, " and Confidence Level = ", format(conf.level, 
                digits = digits), sep = "")
    }, `half.width & n` = {
        if (x.var == "n") {
            x <- seq(min.x, max.x, length = n.points)
            if (sample.type == "two.sample") {
                if (n2.constrained) {
                  if (is.null(xlab)) xlab <- paste("Sample Size (n1), With n2 =", 
                    n2)
                } else {
                  n2 <- x
                  if (is.null(xlab)) xlab <- "Sample Size (n1 and n2)"
                }
            } else if (is.null(xlab)) xlab <- "Sample Size (n)"
            y <- ciNormHalfWidth(n.or.n1 = x, n2 = n2, sigma.hat = sigma.hat, 
                conf.level = conf.level, sample.type = sample.type)
            if (is.null(ylab)) ylab <- "Half-Width"
            if (is.null(main)) main <- paste("Half-Width vs. Sample Size for ", 
                "Confidence Interval for\n", type.string, " with Estimated SD = ", 
                format(sigma.hat, digits = digits), " and Confidence Level = ", 
                format(conf.level, digits = digits), sep = "")
        } else {
            x <- seq(min.x, max.x, length = n.points)
            if (sample.type == "two.sample") {
                if (n2.constrained) {
                  if (is.null(ylab)) ylab <- paste("Sample Size (n1), With n2 =", 
                    n2)
                  y <- ciNormN(half.width = x, sigma.hat = sigma.hat, 
                    conf.level = conf.level, sample.type = "two.sample", 
                    n2 = n2, round.up = round.up, n.max = n.max, 
                    tol = tol, maxiter = maxiter)$n1
                } else {
                  if (is.null(ylab)) ylab <- "Sample Size (n1 and n2)"
                  y <- ciNormN(half.width = x, sigma.hat = sigma.hat, 
                    conf.level = conf.level, sample.type = "two.sample", 
                    round.up = round.up, n.max = n.max, tol = tol, 
                    maxiter = maxiter)
                }
            } else {
                if (is.null(ylab)) ylab <- "Sample Size (n)"
                y <- ciNormN(half.width = x, sigma.hat = sigma.hat, 
                  conf.level = conf.level, sample.type = "one.sample", 
                  round.up = round.up, n.max = n.max, tol = tol, 
                  maxiter = maxiter)
            }
            if (is.null(xlab)) xlab <- "Half-Width"
            if (is.null(main)) main <- paste("Sample Size vs. Half-Width for ", 
                "Confidence Interval for\n", type.string, " with Estimated SD = ", 
                format(sigma.hat, digits = digits), " and Confidence Level = ", 
                format(conf.level, digits = digits), sep = "")
        }
    })
    if (plot.it) {
        if (!add) {
            plot(x, y, type = "n", main = "", sub = "", ..., 
                xlab = xlab, ylab = ylab)
            arg.list <- c(gen.gp.list, list(main = main))
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
