plotTTestDesign <-
function (x.var = "n", y.var = "power", range.x.var = NULL, n.or.n1 = 25, 
    n2 = n.or.n1, delta.over.sigma = switch(alternative, greater = 0.5, 
        less = -0.5, two.sided = ifelse(two.sided.direction == 
            "greater", 0.5, -0.5)), alpha = 0.05, power = 0.95, 
    sample.type = ifelse(!missing(n2), "two.sample", "one.sample"), 
    alternative = "two.sided", two.sided.direction = "greater", 
    approx = FALSE, round.up = FALSE, n.max = 5000, tol = 1e-07, 
    maxiter = 1000, plot.it = TRUE, add = FALSE, n.points = 50, 
    plot.col = "black", plot.lwd = 3 * par("cex"), plot.lty = 1, 
    digits = .Options$digits, ..., main = NULL, xlab = NULL, 
    ylab = NULL, type = "l") 
{
    x.var <- match.arg(x.var, c("n", "delta.over.sigma", "power", 
        "alpha"))
    y.var <- match.arg(y.var, c("power", "delta.over.sigma", 
        "n"))
    sample.type <- match.arg(sample.type, c("one.sample", "two.sample"))
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    two.sided.direction <- match.arg(two.sided.direction, c("greater", 
        "less"))
    if (x.var == y.var) 
        stop("'x.var' and 'y.var' cannot denote the same quantity")
    if (missing(range.x.var)) {
        range.x.var <- switch(x.var, n = c(2, 50), delta.over.sigma = switch(alternative, 
            greater = c(0.5, 2), less = -c(2, 0.5), two.sided = if (two.sided.direction == 
                "greater") c(0.5, 2) else -c(2, 0.5)), power = c(alpha + 
            .Machine$double.eps, 0.95), alpha = c(0.01, 0.2))
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
    if (x.var != "alpha") {
        if (is.null(alpha) || !is.finite(alpha) || !is.vector(alpha, 
            mode = "numeric") || length(alpha) != 1 || alpha <= 
            0 || alpha >= 1) {
            stop("'alpha' must be a scalar greater than 0 and less than 1")
        }
    }
    switch(x.var, n = {
        if (!all(range.x.var == trunc(range.x.var)) || min.x < 
            2) stop(paste("When x.var='n', 'range.x.var' must be", 
            "a vector containing two integers, with the first", 
            "element greater than 1, and the second", "element larger than the first"))
    }, delta.over.sigma = {
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
    else if (sample.type == "two.sample" && n2.constrained && 
        (is.null(n2) || !is.finite(n2) || !is.vector(n2, mode = "numeric") || 
            length(n2) != 1 || n2 < 2 || n2 != trunc(n2))) 
        stop("n2 must be an integer greater than 1")
    if (x.var != "delta.over.sigma" && y.var != "delta.over.sigma" && 
        (is.null(delta.over.sigma) || !is.vector(delta.over.sigma, 
            mode = "numeric") || length(delta.over.sigma) != 
            1)) 
        stop("'delta.over.sigma' must be a numeric scalar")
    if (x.var != "power" && y.var != "power" && (is.null(power) || 
        !is.vector(power, mode = "numeric") || length(power) != 
        1 || power <= 0 || power >= 1)) 
        stop("'power' must be a numeric scalar between 0 and 1")
    if (!is.vector(n.max, mode = "numeric") || length(n.max) != 
        1 || !is.finite(n.max) || n.max != trunc(n.max) || n.max < 
        2) 
        stop("'n.max' must be a positive integer greater than 1")
    if (!is.vector(maxiter, mode = "numeric") || length(maxiter) != 
        1 || !is.finite(maxiter) || maxiter != trunc(maxiter) || 
        maxiter < 2) 
        stop("'maxiter' must be a positive integer greater than 1")
    type.string <- ifelse(sample.type == "one.sample", "One-Sample", 
        "Two-Sample")
    delta.string <- ifelse(sample.type == "one.sample", "(mu - mu0) / sigma", 
        "(mu1 - mu2) / sigma")
    alt.string <- switch(alternative, two.sided = "(Two-Sided Alternative)", 
        less = "(Lower One-Sided Alternative)", greater = "(Upper One-Sided Alternative)")
    if (x.var != "n" && y.var != "n") {
        if (sample.type == "two.sample") 
            n.string <- paste("n1 = ", n.or.n1, ", n2 = ", n2, 
                ",", sep = "")
        else n.string <- paste("n =", n.or.n1)
    }
    if (plot.it) 
        gen.gp.list <- checkGraphicsPars(...)$gen.gp.list
    combo <- paste(sort(c(x.var, y.var)), collapse = " & ")
    switch(combo, `alpha & delta.over.sigma` = {
        x <- seq(min.x, max.x, length = n.points)
        y <- tTestScaledMdd(n.or.n1 = n.or.n1, n2 = n2, alpha = x, 
            power = power, sample.type = sample.type, alternative = alternative, 
            two.sided.direction = two.sided.direction, approx = approx, 
            tol = tol, maxiter = maxiter)
        if (is.null(ylab)) ylab <- delta.string
        if (is.null(xlab)) xlab <- "Alpha"
        if (is.null(main)) main <- paste("Delta/Sigma vs. Alpha for ", 
            type.string, " t-Test with\n", n.string, " and Power = ", 
            format(power, digits = digits), " ", alt.string, 
            sep = "")
    }, `alpha & n` = {
        if (alternative == "greater" && delta.over.sigma <= 0) stop(paste("When alternative='greater',", 
            "'delta.over.sigma' must be positive."))
        if (alternative == "less" && delta.over.sigma >= 0) stop(paste("When alternative='less',", 
            "'delta.over.sigma' must be negative."))
        x <- seq(min.x, max.x, length = n.points)
        if (sample.type == "two.sample") {
            if (n2.constrained) {
                if (is.null(ylab)) ylab <- paste("Sample Size (n1), With n2 =", 
                  n2)
                y <- tTestN(delta.over.sigma = delta.over.sigma, 
                  alpha = x, power = power, sample.type = "two.sample", 
                  alternative = alternative, approx = approx, 
                  n2 = n2, round.up = round.up, n.max = n.max, 
                  tol = tol, maxiter = maxiter)$n1
            } else {
                if (is.null(ylab)) ylab <- "Sample Size (n1 and n2)"
                y <- tTestN(delta.over.sigma = delta.over.sigma, 
                  alpha = x, power = power, sample.type = "two.sample", 
                  alternative = alternative, approx = approx, 
                  round.up = round.up, n.max = n.max, tol = tol, 
                  maxiter = maxiter)
            }
        } else {
            if (is.null(ylab)) ylab <- "Sample Size (n)"
            y <- tTestN(delta.over.sigma = delta.over.sigma, 
                alpha = x, power = power, sample.type = "one.sample", 
                alternative = alternative, approx = approx, round.up = round.up, 
                n.max = n.max, tol = tol, maxiter = maxiter)
        }
        if (is.null(xlab)) xlab <- "Alpha"
        if (is.null(main)) main <- paste("Sample Size vs. Alpha for", 
            type.string, "t-Test with\nDelta/Sigma =", format(delta.over.sigma, 
                digits = digits), "and Power =", format(power, 
                digits = digits), alt.string)
    }, `alpha & power` = {
        x <- seq(min.x, max.x, length = n.points)
        y <- tTestPower(n.or.n1 = n.or.n1, n2 = n2, delta.over.sigma = delta.over.sigma, 
            alpha = x, sample.type = sample.type, alternative = alternative, 
            approx = approx)
        if (is.null(xlab)) xlab <- "Alpha"
        if (is.null(ylab)) ylab <- "Power"
        if (is.null(main)) main <- paste("Power vs. Alpha for ", 
            type.string, " t-Test with\n", n.string, " and Delta/Sigma = ", 
            format(delta.over.sigma, digits = digits), " ", alt.string, 
            sep = "")
    }, `delta.over.sigma & n` = {
        if (x.var == "delta.over.sigma") {
            x <- seq(min.x, max.x, length = n.points)
            if (sample.type == "two.sample") {
                if (n2.constrained) {
                  if (is.null(ylab)) ylab <- paste("Sample Size (n1), With n2 =", 
                    n2)
                  y <- tTestN(delta.over.sigma = x, alpha = alpha, 
                    power = power, sample.type = "two.sample", 
                    alternative = alternative, approx = approx, 
                    n2 = n2, round.up = round.up, n.max = n.max, 
                    tol = tol, maxiter = maxiter)$n1
                } else {
                  if (is.null(ylab)) ylab <- "Sample Size (n1 and n2)"
                  y <- tTestN(delta.over.sigma = x, alpha = alpha, 
                    power = power, sample.type = "two.sample", 
                    alternative = alternative, approx = approx, 
                    round.up = round.up, n.max = n.max, tol = tol, 
                    maxiter = maxiter)
                }
            } else {
                if (is.null(ylab)) ylab <- "Sample Size (n)"
                y <- tTestN(delta.over.sigma = x, alpha = alpha, 
                  power = power, sample.type = "one.sample", 
                  alternative = alternative, approx = approx, 
                  round.up = round.up, n.max = n.max, tol = tol, 
                  maxiter = maxiter)
            }
            if (is.null(xlab)) xlab <- delta.string
            if (is.null(main)) main <- paste("Sample Size vs. Delta/Sigma for", 
                type.string, "t-Test with\nPower =", format(power, 
                  digits = digits), "and Alpha =", format(alpha, 
                  digits = digits), alt.string)
        } else {
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
            y <- tTestScaledMdd(n.or.n1 = x, n2 = n2, alpha = alpha, 
                power = power, sample.type = sample.type, alternative = alternative, 
                two.sided.direction = two.sided.direction, approx = approx, 
                tol = tol, maxiter = maxiter)
            if (is.null(ylab)) ylab <- delta.string
            if (is.null(main)) main <- paste("Delta/Sigma vs. Sample Size for", 
                type.string, "t-Test with\nPower =", format(power, 
                  digits = digits), "and Alpha =", format(alpha, 
                  digits = digits), alt.string)
        }
    }, `delta.over.sigma & power` = {
        if (x.var == "delta.over.sigma") {
            x <- seq(min.x, max.x, length = n.points)
            y <- tTestPower(n.or.n1 = n.or.n1, n2 = n2, delta.over.sigma = x, 
                alpha = alpha, sample.type = sample.type, alternative = alternative, 
                approx = approx)
            if (is.null(xlab)) xlab <- delta.string
            if (is.null(ylab)) ylab <- "Power"
            if (is.null(main)) main <- paste("Power vs. Delta/Sigma for ", 
                type.string, " t-Test with\n", n.string, " and Alpha = ", 
                format(alpha, digits = digits), " ", alt.string, 
                sep = "")
        } else {
            x <- seq(min.x, max.x, length = n.points)
            y <- tTestScaledMdd(n.or.n1 = n.or.n1, n2 = n2, alpha = alpha, 
                power = x, sample.type = sample.type, alternative = alternative, 
                two.sided.direction = two.sided.direction, approx = approx, 
                tol = tol, maxiter = maxiter)
            if (is.null(ylab)) ylab <- delta.string
            if (is.null(xlab)) xlab <- "Power"
            if (is.null(main)) main <- paste("Delta/Sigma vs. Power for ", 
                type.string, " t-Test with\n", n.string, " and Alpha = ", 
                format(alpha, digits = digits), " ", alt.string, 
                sep = "")
        }
    }, `n & power` = {
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
            y <- tTestPower(n.or.n1 = x, n2 = n2, delta.over.sigma = delta.over.sigma, 
                alpha = alpha, sample.type = sample.type, alternative = alternative, 
                approx = approx)
            if (is.null(ylab)) ylab <- "Power"
            if (is.null(main)) main <- paste("Power vs. Sample Size for", 
                type.string, "t-Test with\nDelta/Sigma =", format(delta.over.sigma, 
                  digits = digits), "and Alpha =", format(alpha, 
                  digits = digits), alt.string)
        } else {
            if (alternative == "greater" && delta.over.sigma <= 
                0) stop(paste("When alternative='greater',", 
                "'delta.over.sigma' must be positive."))
            if (alternative == "less" && delta.over.sigma >= 
                0) stop(paste("When alternative='less',", "'delta.over.sigma' must be negative."))
            x <- seq(min.x, max.x, length = n.points)
            if (sample.type == "two.sample") {
                if (n2.constrained) {
                  if (is.null(ylab)) ylab <- paste("Sample Size (n1), With n2 =", 
                    n2)
                  y <- tTestN(delta.over.sigma = delta.over.sigma, 
                    alpha = alpha, power = x, sample.type = "two.sample", 
                    alternative = alternative, approx = approx, 
                    n2 = n2, round.up = round.up, n.max = n.max, 
                    tol = tol, maxiter = maxiter)$n1
                } else {
                  if (is.null(ylab)) ylab <- "Sample Size (n1 and n2)"
                  y <- tTestN(delta.over.sigma = delta.over.sigma, 
                    alpha = alpha, power = x, sample.type = "two.sample", 
                    alternative = alternative, approx = approx, 
                    round.up = round.up, n.max = n.max, tol = tol, 
                    maxiter = maxiter)
                }
            } else {
                if (is.null(ylab)) ylab <- "Sample Size (n)"
                y <- tTestN(delta.over.sigma = delta.over.sigma, 
                  alpha = alpha, power = x, sample.type = "one.sample", 
                  alternative = alternative, approx = approx, 
                  round.up = round.up, n.max = n.max, tol = tol, 
                  maxiter = maxiter)
            }
            if (is.null(xlab)) xlab <- "Power"
            if (is.null(main)) main <- paste("Sample Size vs. Power for", 
                type.string, "t-Test with\nDelta/Sigma =", format(delta.over.sigma, 
                  digits = digits), "and Alpha =", format(alpha, 
                  digits = digits), alt.string)
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
