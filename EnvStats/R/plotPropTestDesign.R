plotPropTestDesign <-
function (x.var = "n", y.var = "power", range.x.var = NULL, n.or.n1 = 25, 
    n2 = n.or.n1, ratio = 1, p.or.p1 = switch(alternative, greater = 0.6, 
        less = 0.4, two.sided = ifelse(two.sided.direction == 
            "greater", 0.6, 0.4)), p0.or.p2 = 0.5, alpha = 0.05, 
    power = 0.95, sample.type = ifelse(!missing(n2) || !missing(ratio), 
        "two.sample", "one.sample"), alternative = "two.sided", 
    two.sided.direction = "greater", approx = TRUE, correct = sample.type == 
        "two.sample", round.up = FALSE, warn = TRUE, n.min = 2, 
    n.max = 10000, tol.alpha = 0.1 * alpha, tol = 1e-07, maxiter = 1000, 
    plot.it = TRUE, add = FALSE, n.points = 50, plot.col = "black", 
    plot.lwd = 3 * par("cex"), plot.lty = 1, digits = .Options$digits, 
    cex.main = par("cex"), ..., main = NULL, xlab = NULL, ylab = NULL, 
    type = "l") 
{
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    two.sided.direction <- match.arg(two.sided.direction, c("greater", 
        "less"))
    sample.type <- match.arg(sample.type, c("one.sample", "two.sample"))
    if (sample.type == "two.sample" && !approx) 
        stop(paste("Exact method not available for", "the two-sample case.  Set approx=T."))
    x.var <- match.arg(x.var, c("n", "delta", "power", "alpha"))
    y.var <- match.arg(y.var, c("power", "delta", "n"))
    if (x.var == y.var) 
        stop("'x.var' and 'y.var' cannot denote the same quantity")
    if (missing(range.x.var)) {
        range.x.var <- switch(x.var, n = c(20, 400), delta = switch(alternative, 
            greater = c(0.05, 0.2), less = -c(0.2, 0.05), two.sided = if (two.sided.direction == 
                "greater") c(0.05, 0.2) else -c(0.2, 0.05)), 
            power = c(alpha + .Machine$double.eps, 0.95), alpha = c(0.01, 
                0.2))
    }
    else {
        if (is.null(range.x.var) || !all(is.finite(range.x.var)) || 
            !is.vector(range.x.var, mode = "numeric") || length(range.x.var) != 
            2 || range.x.var[1] >= range.x.var[2]) 
            stop(paste("'range.x.var' must be a numeric vector of length 2", 
                "with no missing (NA), infinite(-Inf, Inf), or undefined(NaN) values", 
                "and the first element must be less than the second"))
    }
    min.x <- range.x.var[1]
    max.x <- range.x.var[2]
    if (!is.vector(n.points, mode = "numeric") || length(n.points) != 
        1 || n.points != trunc(n.points) || n.points < 2) 
        stop("'n.points' must be an integer larger than 1")
    if (x.var != "alpha") {
        if (is.null(alpha) || !is.finite(alpha) || !is.vector(alpha, 
            mode = "numeric") || length(alpha) != 1 || alpha <= 
            0 || alpha >= 1) {
            stop("'alpha' must be a scalar between 0 and 1")
        }
    }
    switch(x.var, n = {
        if (!all(range.x.var == trunc(range.x.var)) || min.x < 
            2) stop(paste("When x.var='n', 'range.x.var' must be", 
            "a vector containing two integers, with the first", 
            "element greater than 1, and the second", "element larger than the first"))
    }, delta = {
        if (alternative == "greater" && (min.x <= 0 || min.x >= 
            1)) stop(paste("When alternative='greater',", "and x.var='delta', all values of", 
            "'delta' must be positive and less than 1,", "so the first element of 'range.x.var' must be", 
            "positive and less than 1."))
        if (alternative == "less" && (max.x >= 0 || max.x <= 
            -1)) stop(paste("When alternative='less',", "and x.var='delta', all values of", 
            "'delta' must be negative and greater than -1,", 
            "so the second element of 'range.x.var' must be", 
            "negative and greater than 1."))
    }, power = {
        if (min.x < alpha || max.x >= 1) stop(paste("When x.var='power', 'range.x.var' must be", 
            "a vector of two elements, with the first element", 
            "greater than or equal to 'alpha', the second element", 
            "less than 1, and the second", "element larger than the first."))
    }, alpha = {
        if (min.x <= 0 || max.x >= 1) stop(paste("When x.var='alpha', 'range.x.var' must be", 
            "a vector of two elements, with the first element", 
            "greater than 0, the second element less than 1, and the second", 
            "element larger than the first."))
    })
    n2.constrained <- !missing(n2) && !is.null(n2)
    ratio.supplied <- !missing(ratio) && !is.null(ratio)
    if (x.var != "n" && y.var != "n") {
        if (is.null(n.or.n1) || !is.finite(n.or.n1) || !is.vector(n.or.n1, 
            mode = "numeric") || length(n.or.n1) != 1 || n.or.n1 < 
            2 || n.or.n1 != trunc(n.or.n1)) 
            stop("'n.or.n1 must be an integer greater than 1")
        if (sample.type == "two.sample" && (is.null(n2) || !is.finite(n2) || 
            !is.vector(n2, mode = "numeric") || length(n2) != 
            1 || n2 < 2 || n2 != trunc(n2))) 
            stop("n2 must be an integer greater than 1")
    }
    else {
        if (sample.type == "two.sample") {
            if (x.var == "n") {
                if (n2.constrained && (is.null(n2) || !is.finite(n2) || 
                  !is.vector(n2, mode = "numeric") || length(n2) != 
                  1 || n2 < 2 || n2 != trunc(n2))) 
                  stop("n2 must be an integer greater than 1")
            }
            else {
                if (ratio.supplied) {
                  if (!is.vector(ratio, mode = "numeric") || 
                    length(ratio) != 1 || ratio < 1) 
                    stop("'ratio' must be a numeric scalar greater than or equal to 1")
                  ratio.gt.1 <- ratio > 1
                }
                else ratio.gt.1 <- FALSE
            }
        }
    }
    if (!is.vector(p0.or.p2, mode = "numeric") || length(p0.or.p2) != 
        1 || p0.or.p2 <= 0 || p0.or.p2 >= 1) 
        stop("'p0.or.p2' must be a numeric scalar greater than 0 and less than 1")
    if (x.var != "delta" && y.var != "delta") {
        if (!is.vector(p.or.p1, mode = "numeric") || length(p.or.p1) != 
            1 || p.or.p1 <= 0 || p.or.p1 >= 1) 
            stop("'p.or.p1' must be a numeric scalar greater than 0 and less than 1")
        if (p.or.p1 == p0.or.p2) 
            stop("'p.or.p1' cannot equal 'p0.or.p2'")
    }
    if (x.var != "power" && y.var != "power") {
        if (is.null(power) || !is.vector(power, mode = "numeric") || 
            length(power) != 1 || power <= 0 || power >= 1) 
            stop("'power' must be a numeric scalar greater than 0 and less than 1")
        if (x.var != "alpha" && power <= alpha) 
            stop("power must be greater than alpha")
    }
    type.string <- paste(ifelse(sample.type == "one.sample", 
        "One-Sample", "Two-Sample"), "Proportion Test with")
    delta.string <- ifelse(sample.type == "one.sample", "(p - p0)", 
        "(p1 - p2)")
    p0.or.p2.string <- paste(ifelse(sample.type == "one.sample", 
        "p0", "p2"), "=", format(p0.or.p2, digits = digits))
    p.or.p1.string <- paste(ifelse(sample.type == "one.sample", 
        "p", "p1"), "=", format(p.or.p1, digits = digits))
    alt.string <- switch(alternative, two.sided = "(Two-Sided Alternative)", 
        less = "(Lower One-Sided Alternative)", greater = "(Upper One-Sided Alternative)")
    if (sample.type == "two.sample") 
        n.string <- paste("n1 = ", n.or.n1, ", n2 = ", n2, sep = "")
    else n.string <- paste("n =", n.or.n1)
    alpha.string <- paste(ifelse(approx, "Alpha =", "Alpha <="), 
        format(alpha, digits = digits))
    power.string <- paste("Power =", format(power, digits = digits))
    approx.string <- paste("Based on", ifelse(approx, "Normal Approximation", 
        "Exact Test"))
    if (plot.it) 
        gen.gp.list <- checkGraphicsPars(...)$gen.gp.list
    combo <- paste(sort(c(x.var, y.var)), collapse = " & ")
    switch(combo, `alpha & delta` = {
        x <- seq(range.x.var[1], range.x.var[2], length = n.points)
        y <- propTestMdd(n.or.n1 = n.or.n1, n2 = n2, p0.or.p2 = p0.or.p2, 
            alpha = x, power = power, sample.type = sample.type, 
            alternative = alternative, two.sided.direction = two.sided.direction, 
            approx = approx, correct = correct, tol = tol, warn = warn, 
            return.exact.list = FALSE)
        if (is.null(ylab)) ylab <- delta.string
        if (is.null(xlab)) xlab <- "Alpha"
        line1 <- paste("Delta vs. Alpha for", type.string)
        line2 <- paste(n.string, ", ", p0.or.p2.string, ", and ", 
            power.string, " ", alt.string, sep = "")
    }, `alpha & n` = {
        if (alternative == "greater" && p.or.p1 <= p0.or.p2) stop(paste("When alternative='greater',", 
            "'p.or.p1' must be greater than 'p0.or.p2'"))
        if (alternative == "less" && p.or.p1 >= p0.or.p2) stop(paste("When alternative='less',", 
            "'p.or.p1' must be less than 'p0.or.p2'"))
        x <- seq(range.x.var[1], range.x.var[2], length = n.points)
        if (sample.type == "two.sample") {
            if (ratio.gt.1) {
                if (is.null(ylab)) ylab <- paste("Sample Size (n1), With n2 =", 
                  ratio, "* n1")
                y <- propTestN(p.or.p1 = p.or.p1, p0.or.p2 = p0.or.p2, 
                  alpha = x, power = power, sample.type = "two.sample", 
                  alternative = alternative, ratio = ratio, approx = approx, 
                  correct = correct, round.up = round.up, warn = warn, 
                  return.exact.list = FALSE)$n1
            } else {
                if (is.null(ylab)) ylab <- "Sample Size (n1 and n2)"
                y <- propTestN(p.or.p1 = p.or.p1, p0.or.p2 = p0.or.p2, 
                  alpha = x, power = power, sample.type = "two.sample", 
                  alternative = alternative, ratio = 1, approx = approx, 
                  correct = correct, round.up = round.up, warn = warn, 
                  return.exact.list = FALSE)
            }
        } else {
            if (is.null(ylab)) ylab <- "Sample Size (n)"
            y <- propTestN(p.or.p1 = p.or.p1, p0.or.p2 = p0.or.p2, 
                alpha = x, power = power, sample.type = "one.sample", 
                alternative = alternative, approx = approx, correct = correct, 
                round.up = round.up, warn = warn, return.exact.list = FALSE, 
                tol = tol, n.min = n.min, n.max = n.max, tol.alpha = tol.alpha, 
                maxiter = maxiter)
        }
        if (is.null(xlab)) xlab <- "Alpha"
        line1 <- paste("Sample Size vs. Alpha for", type.string)
        line2 <- paste(p.or.p1.string, ", ", p0.or.p2.string, 
            ", and ", power.string, " ", alt.string, sep = "")
    }, `alpha & power` = {
        x <- seq(range.x.var[1], range.x.var[2], length = n.points)
        y <- propTestPower(n.or.n1 = n.or.n1, p.or.p1 = p.or.p1, 
            n2 = n2, p0.or.p2 = p0.or.p2, alpha = x, sample.type = sample.type, 
            alternative = alternative, approx = approx, correct = correct, 
            warn = warn, return.exact.list = FALSE)
        if (is.null(xlab)) xlab <- "Alpha"
        if (is.null(ylab)) ylab <- "Power"
        line1 <- paste("Power vs. Alpha for", type.string)
        line2 <- paste(n.string, ", ", p.or.p1.string, ", and ", 
            p0.or.p2.string, " ", alt.string, sep = "")
    }, `delta & n` = {
        if (x.var == "delta") {
            if (any(range.x.var + p0.or.p2 <= 0) || any(range.x.var + 
                p0.or.p2 >= 1)) stop(paste("When x.var='delta' and y.var='n',", 
                "both elements of", "'range.x.var' + 'p0.or.p2' must fall in the", 
                "interval (0,1)"))
            x <- seq(range.x.var[1], range.x.var[2], length = n.points)
            if (sample.type == "two.sample") {
                if (ratio.gt.1) {
                  if (is.null(ylab)) ylab <- paste("Sample Size (n1), With n2 =", 
                    ratio, "* n1")
                  y <- propTestN(p.or.p1 = x + p0.or.p2, p0.or.p2 = p0.or.p2, 
                    alpha = alpha, power = power, sample.type = "two.sample", 
                    alternative = alternative, ratio = ratio, 
                    approx = approx, correct = correct, round.up = round.up, 
                    warn = warn, return.exact.list = FALSE)$n1
                } else {
                  if (is.null(ylab)) ylab <- "Sample Size (n1 and n2)"
                  y <- propTestN(p.or.p1 = x + p0.or.p2, p0.or.p2 = p0.or.p2, 
                    alpha = alpha, power = power, sample.type = "two.sample", 
                    alternative = alternative, approx = approx, 
                    correct = correct, round.up = round.up, warn = warn, 
                    return.exact.list = FALSE)
                }
            } else {
                if (is.null(ylab)) ylab <- "Sample Size (n)"
                y <- propTestN(p.or.p1 = x + p0.or.p2, p0.or.p2 = p0.or.p2, 
                  alpha = alpha, power = power, sample.type = "one.sample", 
                  alternative = alternative, approx = approx, 
                  correct = correct, round.up = round.up, warn = warn, 
                  return.exact.list = FALSE, tol = tol, n.min = n.min, 
                  n.max = n.max, tol.alpha = tol.alpha, maxiter = maxiter)
            }
            if (is.null(xlab)) xlab <- delta.string
            line1 <- paste("Sample Size vs. Delta for", type.string)
            line2 <- paste(p0.or.p2.string, ", ", power.string, 
                ", and ", alpha.string, " ", alt.string, sep = "")
        } else {
            if (approx) x <- seq(range.x.var[1], range.x.var[2], 
                length = n.points) else x <- seq(min.x, max.x, 
                by = ceiling((max.x - min.x + 1)/n.points))
            if (sample.type == "two.sample") {
                if (n2.constrained) {
                  if (is.null(xlab)) xlab <- paste("Sample Size (n1), With n2 =", 
                    n2)
                } else {
                  n2 <- x
                  if (is.null(xlab)) xlab <- "Sample Size (n1 and n2)"
                }
            } else if (is.null(xlab)) xlab <- "Sample Size (n)"
            y <- propTestMdd(n.or.n1 = x, n2 = n2, p0.or.p2 = p0.or.p2, 
                alpha = alpha, power = power, sample.type = sample.type, 
                alternative = alternative, two.sided.direction = two.sided.direction, 
                approx = approx, correct = correct, tol = tol, 
                warn = warn, return.exact.list = FALSE)
            if (is.null(ylab)) ylab <- delta.string
            line1 <- paste("Delta vs. Sample Size for", type.string)
            line2 <- paste(p0.or.p2.string, ", ", power.string, 
                ", and ", alpha.string, " ", alt.string, sep = "")
        }
    }, `delta & power` = {
        if (x.var == "delta") {
            if (any(range.x.var + p0.or.p2 <= 0) || any(range.x.var + 
                p0.or.p2 >= 1)) stop(paste("When x.var='delta' and y.var='power',", 
                "both elements of", "'range.x.var' + 'p0.or.p2' must fall in the", 
                "interval (0,1)"))
            x <- seq(range.x.var[1], range.x.var[2], length = n.points)
            y <- propTestPower(n.or.n1 = n.or.n1, p.or.p1 = x + 
                p0.or.p2, n2 = n2, p0.or.p2 = p0.or.p2, alpha = alpha, 
                sample.type = sample.type, alternative = alternative, 
                approx = approx, correct = correct, warn = warn, 
                return.exact.list = FALSE)
            if (is.null(xlab)) xlab <- delta.string
            if (is.null(ylab)) ylab <- "Power"
            line1 <- paste("Power vs. Delta for", type.string)
            line2 <- paste(n.string, ", ", p0.or.p2.string, ", and ", 
                alpha.string, " ", alt.string, sep = "")
        } else {
            x <- seq(range.x.var[1], range.x.var[2], length = n.points)
            y <- propTestMdd(n.or.n1 = n.or.n1, n2 = n2, p0.or.p2 = p0.or.p2, 
                alpha = alpha, power = x, sample.type = sample.type, 
                alternative = alternative, two.sided.direction = two.sided.direction, 
                approx = approx, correct = correct, tol = tol, 
                warn = warn, return.exact.list = FALSE)
            if (is.null(ylab)) ylab <- delta.string
            if (is.null(xlab)) xlab <- "Power"
            line1 <- paste("Delta vs. Power for", type.string)
            line2 <- paste(n.string, ", ", p0.or.p2.string, ", and ", 
                alpha.string, " ", alt.string, sep = "")
        }
    }, `n & power` = {
        if (x.var == "n") {
            if (approx) x <- seq(range.x.var[1], range.x.var[2], 
                length = n.points) else x <- seq(min.x, max.x, 
                by = ceiling((max.x - min.x + 1)/n.points))
            if (sample.type == "two.sample") {
                if (n2.constrained) {
                  if (is.null(xlab)) xlab <- paste("Sample Size (n1), With n2 =", 
                    n2)
                } else {
                  n2 <- x
                  if (is.null(xlab)) xlab <- "Sample Size (n1 and n2)"
                }
            } else if (is.null(xlab)) xlab <- "Sample Size (n)"
            y <- propTestPower(n.or.n1 = x, p.or.p1 = p.or.p1, 
                n2 = n2, p0.or.p2 = p0.or.p2, alpha = alpha, 
                sample.type = sample.type, alternative = alternative, 
                approx = approx, correct = correct, warn = warn, 
                return.exact.list = FALSE)
            if (is.null(ylab)) ylab <- "Power"
            line1 <- paste("Power vs. Sample Size for", type.string)
            line2 <- paste(p.or.p1.string, ", ", p0.or.p2.string, 
                ", and ", alpha.string, " ", alt.string, sep = "")
        } else {
            if (alternative == "greater" && p.or.p1 <= p0.or.p2) stop(paste("When alternative='greater',", 
                "'p.or.p1' must be greater than 'p0.or.p2'"))
            if (alternative == "less" && p.or.p1 >= p0.or.p2) stop(paste("When alternative='less',", 
                "'p.or.p1' must be less than 'p0.or.p2'"))
            x <- seq(range.x.var[1], range.x.var[2], length = n.points)
            if (sample.type == "two.sample") {
                if (ratio.gt.1) {
                  if (is.null(ylab)) ylab <- paste("Sample Size (n1), With n2 =", 
                    ratio, "* n1")
                  y <- propTestN(p.or.p1 = p.or.p1, p0.or.p2 = p0.or.p2, 
                    alpha = alpha, power = x, sample.type = "two.sample", 
                    alternative = alternative, ratio = ratio, 
                    approx = approx, correct = correct, round.up = round.up, 
                    warn = warn, return.exact.list = FALSE)$n1
                } else {
                  if (is.null(ylab)) ylab <- "Sample Size (n1 and n2)"
                  y <- propTestN(p.or.p1 = p.or.p1, p0.or.p2 = p0.or.p2, 
                    alpha = alpha, power = x, sample.type = "two.sample", 
                    alternative = alternative, approx = approx, 
                    correct = correct, round.up = round.up, warn = warn, 
                    return.exact.list = FALSE)
                }
            } else {
                if (is.null(ylab)) ylab <- "Sample Size (n)"
                y <- propTestN(p.or.p1 = p.or.p1, p0.or.p2 = p0.or.p2, 
                  alpha = alpha, power = x, sample.type = "one.sample", 
                  alternative = alternative, approx = approx, 
                  correct = correct, round.up = round.up, warn = warn, 
                  return.exact.list = FALSE, tol = tol, n.min = n.min, 
                  n.max = n.max, tol.alpha = tol.alpha, maxiter = maxiter)
            }
            if (is.null(xlab)) xlab <- "Power"
            line1 <- paste("Sample Size vs. Power for", type.string)
            line2 <- paste(p.or.p1.string, ", ", p0.or.p2.string, 
                ", and ", alpha.string, " ", alt.string, sep = "")
        }
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
                mtext(text = approx.string, side = 3, line = 0.5, 
                  cex = cex.main)
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
