plotCiBinomDesign <-
function (x.var = "n", y.var = "half.width", range.x.var = NULL, 
    n.or.n1 = 25, p.hat.or.p1.hat = 0.5, n2 = n.or.n1, p2.hat = 0.4, 
    ratio = 1, half.width = 0.05, conf.level = 0.95, sample.type = "one.sample", 
    ci.method = "score", correct = TRUE, warn = TRUE, n.or.n1.min = 2, 
    n.or.n1.max = 10000, tol.half.width = 0.005, tol.p.hat = 0.005, 
    maxiter = 10000, plot.it = TRUE, add = FALSE, n.points = 100, 
    plot.col = 1, plot.lwd = 3 * par("cex"), plot.lty = 1, digits = .Options$digits, 
    main = NULL, xlab = NULL, ylab = NULL, type = "l", ...) 
{
    x.var <- match.arg(x.var, c("n", "half.width", "p.hat", "conf.level"))
    y.var <- match.arg(y.var, c("half.width", "n"))
    if (x.var == y.var) 
        stop("'x.var' and 'y.var' cannot denote the same quantity")
    if (missing(range.x.var)) {
        range.x.var <- switch(x.var, n = c(10, 50), half.width = c(0.03, 
            0.1), p.hat = c(0.5, 0.9), conf.level = c(0.8, 0.99))
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
    if (missing(sample.type)) {
        sample.type <- ifelse(!missing(n2) || !missing(p2.hat) || 
            !missing(ratio), "two.sample", "one.sample")
    }
    sample.type <- match.arg(sample.type, c("one.sample", "two.sample"))
    ci.method <- match.arg(ci.method, c("score", "exact", "adjusted Wald", 
        "Wald"))
    switch(x.var, n = {
        if (!all(range.x.var == trunc(range.x.var)) || min.x < 
            2) stop(paste("When x.var='n', 'range.x.var' must be", 
            "a vector containing two integers, with the first", 
            "element greater than 1, and the second", "element larger than the first."))
    }, half.width = {
        if (sample.type == "one.sample" && (range.x.var[1] < 
            .Machine$double.eps || range.x.var[2] >= 0.5 || range.x.var[1] >= 
            range.x.var[2])) stop(paste("When x.var='half.width' and sample.type='one.sample',", 
            "'range.x.var' must be", "a vector of two elements, with the first", 
            "element greater than 0, the second element less than 0.5,", 
            "and the second element larger than the first."))
        if (sample.type == "two.sample" && (range.x.var[1] < 
            .Machine$double.eps || range.x.var[2] > 1 - .Machine$double.eps || 
            range.x.var[1] >= range.x.var[2])) stop(paste("When x.var='half.width' and sample.type='one.sample',", 
            "'range.x.var' must be", "a vector of two elements, with the first", 
            "element greater than 0, the second element less than 1,", 
            "and the second element larger than the first."))
    }, p.hat = {
        if (range.x.var[1] < .Machine$double.eps || range.x.var[2] > 
            1 - .Machine$double.eps || range.x.var[1] >= range.x.var[2]) stop(paste("When x.var='p.hat', 'range.x.var' must be", 
            "a vector containing two elements, with the first element", 
            "greater than 0, the second element less than 1,", 
            "and with the second element larger than the first"))
    }, conf.level = {
        if (range.x.var[1] < .Machine$double.eps || range.x.var[2] > 
            1 - .Machine$double.eps || range.x.var[1] >= range.x.var[2]) stop(paste("When x.var='conf.level', 'range.x.var' must be", 
            "a vector of two elements, with the first element", 
            "greater than 0, the second element less than 1, and the second", 
            "element larger than the first."))
    })
    if (x.var != "n" && y.var != "n") {
        if (!is.vector(n.or.n1, mode = "numeric") || length(n.or.n1) != 
            1 || n.or.n1 != trunc(n.or.n1) || n.or.n1 < 2) 
            stop("'n.or.n1' must be an integer greater than 1")
    }
    if (x.var != "half.width" && y.var != "half.width") {
        if (!is.vector(half.width, mode = "numeric") || length(half.width) != 
            1 || half.width < .Machine$double.eps) 
            stop("'half.width' must be a positive numeric scalar")
    }
    if (x.var != "p.hat") {
        if (!is.vector(p.hat.or.p1.hat, mode = "numeric") || 
            length(p.hat.or.p1.hat) != 1 || p.hat.or.p1.hat < 
            .Machine$double.eps || p.hat.or.p1.hat > 1 - .Machine$double.eps) 
            stop("'p.hat.or.p1.hat' must be a positive numeric scalar between 0 and 1")
    }
    if (x.var != "conf.level") {
        if (!is.vector(conf.level, mode = "numeric") || length(conf.level) != 
            1 || conf.level < .Machine$double.eps || conf.level > 
            1 - .Machine$double.eps) 
            stop("'conf.level' must be a numeric scalar between 0 and 1")
    }
    if (sample.type == "two.sample") {
        if (ci.method == "exact") 
            stop("Exact method not available for two-sample confidence interval")
        if (x.var == "n" || y.var == "n") {
            if (!is.vector(ratio, mode = "numeric") || length(ratio) != 
                1 || any(ratio < 1)) 
                stop("'ratio' must be a numeric scalar greater than or equal to 1")
        }
        else if (!is.vector(n2, mode = "numeric") || length(n2) != 
            1 || n2 != trunc(n2) || n2 < 2) {
            stop("'n2' must be an integer greater than 1")
        }
        if (!is.vector(p2.hat, mode = "numeric") || length(p2.hat) != 
            1 || p2.hat < .Machine$double.eps || p2.hat > 1 - 
            .Machine$double.eps) 
            stop("'p2.hat' must be a positive numeric scalar between 0 and 1")
    }
    if (y.var == "n") {
        if (length(n.or.n1.min) != 1 || length(n.or.n1.max) != 
            1 || n.or.n1.min < 2 || n.or.n1.max <= n.or.n1.min || 
            n.or.n1.min != round(n.or.n1.min) || n.or.n1.max != 
            round(n.or.n1.max)) 
            stop(paste("'n.or.n1.min' must be an integer greater than 1", 
                "and 'n.or.n1.max' must be an integer greater than 'n.or.n1.min'"))
        if (length(tol.half.width) != 1 || tol.half.width <= 
            .Machine$double.eps || tol.half.width >= 0.5) 
            stop("'tol.half.width' must be a scalar greater than 0 and less than 0.5")
        if (length(tol.p.hat) != 1 || tol.p.hat <= .Machine$double.eps || 
            tol.p.hat >= 0.5) 
            stop("'tol.p.hat' must be a scalar greater than 0 and less than 0.5")
    }
    if (!is.vector(n.points, mode = "numeric") || length(n.points) != 
        1 || n.points != trunc(n.points) || n.points < 2) 
        stop("'n.points' must be an integer larger than 1")
    type.string <- ifelse(sample.type == "one.sample", "p,", 
        "(p1 - p2),")
    if (plot.it) 
        gen.gp.list <- checkGraphicsPars(...)$gen.gp.list
    combo <- paste(sort(c(x.var, y.var)), collapse = " & ")
    p.string <- ifelse(sample.type == "two.sample", paste(", Estimated p1 = ", 
        format(p.hat.or.p1.hat, digits = digits), ", and Estimated p2 = ", 
        format(p2.hat, digits = digits), sep = ""), paste(", and Estimated p =", 
        format(p.hat.or.p1.hat, digits = digits)))
    ratio.gt.1 <- ratio > 1
    n2.constrained <- !missing(n2) && !is.null(n2)
    is.null.main <- is.null(main)
    switch(combo, `conf.level & half.width` = {
        x <- seq(range.x.var[1], range.x.var[2], length = n.points)
        y <- ciBinomHalfWidth(n.or.n1 = n.or.n1, n2 = n2, p.hat.or.p1.hat = p.hat.or.p1.hat, 
            p2.hat = p2.hat, conf.level = x, sample.type = sample.type, 
            ci.method = ci.method, correct = correct, warn = warn)$half.width
        if (is.null(ylab)) ylab <- "Half-Width"
        if (is.null(xlab)) xlab <- "Confidence Level"
        n.string <- ifelse(sample.type == "two.sample", paste("n1 = ", 
            n.or.n1, ", n2 = ", n2, sep = ""), paste("n =", n.or.n1))
        if (is.null(main)) main <- paste("Half-Width vs. Confidence Level for ", 
            "Confidence Interval for ", type.string, "\nwith ", 
            n.string, p.string, sep = "")
    }, `conf.level & n` = {
        x <- seq(range.x.var[1], range.x.var[2], length = n.points)
        if (sample.type == "two.sample") {
            if (is.null(ylab)) {
                ylab <- ifelse(ratio.gt.1, paste("Sample Size (n1), With n2 =", 
                  ratio, "* n1"), "Sample Size (n1 and n2)")
            }
            y <- ciBinomN(half.width = half.width, p.hat.or.p1.hat = p.hat.or.p1.hat, 
                p2.hat = p2.hat, conf.level = x, sample.type = "two.sample", 
                ratio = ratio, ci.method = ci.method, correct = correct, 
                warn = warn, n.or.n1.min = n.or.n1.min, n.or.n1.max = n.or.n1.max, 
                tol.half.width = tol.half.width, tol.p.hat = tol.p.hat, 
                maxiter = maxiter)$n1
        } else {
            if (is.null(ylab)) ylab <- "Sample Size (n)"
            y <- ciBinomN(half.width = half.width, p.hat.or.p1.hat = p.hat.or.p1.hat, 
                conf.level = x, sample.type = "one.sample", ci.method = ci.method, 
                correct = correct, warn = warn, n.or.n1.min = n.or.n1.min, 
                n.or.n1.max = n.or.n1.max, tol.half.width = tol.half.width, 
                tol.p.hat = tol.p.hat, maxiter = maxiter)$n
        }
        if (is.null(xlab)) xlab <- "Confidence Level"
        if (is.null(main)) main <- paste("Sample Size vs. Confidence Level for ", 
            "Confidence Interval for ", type.string, "\nwith Half-Width = ", 
            format(half.width, digits = digits), p.string, sep = "")
    }, `n & p.hat` = {
        x <- seq(range.x.var[1], range.x.var[2], length = n.points)
        if (sample.type == "two.sample") {
            if (is.null(ylab)) {
                ylab <- ifelse(ratio.gt.1, paste("Sample Size (n1), With n2 =", 
                  ratio, "* n1"), "Sample Size (n1 and n2)")
            }
            y <- ciBinomN(half.width = half.width, p.hat.or.p1.hat = x, 
                p2.hat = p2.hat, conf.level = conf.level, sample.type = "two.sample", 
                ratio = ratio, ci.method = ci.method, correct = correct, 
                warn = warn, n.or.n1.min = n.or.n1.min, n.or.n1.max = n.or.n1.max, 
                tol.half.width = tol.half.width, tol.p.hat = tol.p.hat, 
                maxiter = maxiter)$n1
            if (is.null(xlab)) xlab <- "Estimated p1"
            if (is.null(main)) main <- paste("Sample Size vs. Estimated p1 for ", 
                "Confidence Interval for ", type.string, "\nwith Half-Width = ", 
                format(half.width, digits = digits), " Confidence Level = ", 
                format(conf.level, digits = digits), ", and Estimated p2 = ", 
                p2.hat, sep = "")
        } else {
            if (is.null(ylab)) ylab <- "Sample Size (n)"
            y <- ciBinomN(half.width = half.width, p.hat.or.p1.hat = x, 
                conf.level = conf.level, sample.type = "one.sample", 
                ci.method = ci.method, correct = correct, warn = warn, 
                n.or.n1.min = n.or.n1.min, n.or.n1.max = n.or.n1.max, 
                tol.half.width = tol.half.width, tol.p.hat = tol.p.hat, 
                maxiter = maxiter)$n
            if (is.null(xlab)) xlab <- "Estiamted p"
            if (is.null(main)) main <- paste("Sample Size vs. Estimated p for ", 
                "Confidence Interval for ", type.string, "\nwith Half-Width = ", 
                format(half.width, digits = digits), " and Confidence Level = ", 
                format(conf.level, digits = digits), sep = "")
        }
    }, `half.width & p.hat` = {
        x <- seq(range.x.var[1], range.x.var[2], length = n.points)
        if (sample.type == "two.sample") {
            y <- ciBinomHalfWidth(n.or.n1 = n.or.n1, n2 = n2, 
                p.hat.or.p1.hat = x, p2.hat = p2.hat, conf.level = conf.level, 
                sample.type = "two.sample", ci.method = ci.method, 
                correct = correct, warn = warn)$half.width
            if (is.null(xlab)) xlab <- "Estimated p1"
            if (is.null(main)) main <- paste("Half-Width vs. Estimated p1 for ", 
                "Confidence Interval for ", type.string, "\nwith ", 
                "n1 = ", n.or.n1, ", n2 = ", n2, ",", " Estimated p2 = ", 
                p2.hat, ", and Confidence Level = ", format(conf.level, 
                  digits = digits), sep = "")
        } else {
            y <- ciBinomHalfWidth(n.or.n1 = n.or.n1, p.hat.or.p1.hat = x, 
                conf.level = conf.level, sample.type = "one.sample", 
                ci.method = ci.method, correct = correct, warn = warn)$half.width
            if (is.null(xlab)) xlab <- "Estimated p"
            if (is.null(main)) main <- paste("Half-Width vs. Estimated p for ", 
                "Confidence Interval for ", type.string, "\nwith n = ", 
                n.or.n1, " and Confidence Level = ", format(conf.level, 
                  digits = digits), sep = "")
        }
        if (is.null(ylab)) ylab <- "Half-Width"
    }, `half.width & n` = {
        if (x.var == "n") {
            x <- range.x.var[1]:range.x.var[2]
            if (sample.type == "two.sample") {
                if (n2.constrained) {
                  if (is.null(xlab)) xlab <- paste("Sample Size (n1), With n2 =", 
                    n2)
                } else {
                  n2 <- x
                  if (is.null(xlab)) xlab <- "Sample Size (n1 and n2)"
                }
            } else {
                if (is.null(xlab)) xlab <- "Sample Size (n)"
            }
            y <- ciBinomHalfWidth(n.or.n1 = x, n2 = n2, p.hat.or.p1.hat = p.hat.or.p1.hat, 
                p2.hat = p2.hat, conf.level = conf.level, sample.type = sample.type, 
                ci.method = ci.method, correct = correct, warn = warn)$half.width
            if (is.null(ylab)) ylab <- "Half-Width"
            if (is.null(main)) main <- paste("Half-Width vs. Sample Size for ", 
                "Confidence Interval for ", type.string, "\nwith Confidence Level = ", 
                format(conf.level, digits = digits), p.string, 
                sep = "")
        } else {
            x <- seq(range.x.var[1], range.x.var[2], length = n.points)
            if (sample.type == "two.sample") {
                if (is.null(ylab)) {
                  ylab <- ifelse(ratio.gt.1, paste("Sample Size (n1), With n2 =", 
                    ratio, "* n1"), "Sample Size (n1 and n2)")
                }
                y <- ciBinomN(half.width = x, p.hat.or.p1.hat = p.hat.or.p1.hat, 
                  p2.hat = p2.hat, conf.level = conf.level, sample.type = "two.sample", 
                  ratio = ratio, ci.method = ci.method, correct = correct, 
                  warn = warn, n.or.n1.min = n.or.n1.min, n.or.n1.max = n.or.n1.max, 
                  tol.half.width = tol.half.width, tol.p.hat = tol.p.hat, 
                  maxiter = maxiter)$n1
            } else {
                if (is.null(ylab)) ylab <- "Sample Size (n)"
                y <- ciBinomN(half.width = x, p.hat.or.p1.hat = p.hat.or.p1.hat, 
                  conf.level = conf.level, sample.type = "one.sample", 
                  ci.method = ci.method, correct = correct, warn = warn, 
                  n.or.n1.min = n.or.n1.min, n.or.n1.max = n.or.n1.max, 
                  tol.half.width = tol.half.width, tol.p.hat = tol.p.hat, 
                  maxiter = maxiter)$n
            }
            if (is.null(xlab)) xlab <- "Half-Width"
            if (is.null(main)) main <- paste("Sample Size vs. Half-Width for ", 
                "Confidence Interval for ", type.string, "\nwith Confidence Level = ", 
                format(conf.level, digits = digits), p.string, 
                sep = "")
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
            if (is.null.main) {
                method <- switch(ci.method, score = {
                  ifelse(!correct, "Score normal approximation", 
                    "Score normal approximation, with continuity correction")
                }, exact = "Exact", Wald = {
                  ifelse(!correct, "Wald normal approximation", 
                    "Wald normal approximation, with continuity correction")
                }, `adjusted Wald` = "Adjusted Wald normal approximation")
                mtext(paste("CI Method =", method), line = 0, 
                  cex = par("cex"))
            }
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
