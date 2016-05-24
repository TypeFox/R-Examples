plotLinearTrendTestDesign <-
function (x.var = "n", y.var = "power", range.x.var = NULL, n = 12, 
    slope.over.sigma = switch(alternative, greater = 0.1, less = -0.1, 
        two.sided = ifelse(two.sided.direction == "greater", 
            0.1, -0.1)), alpha = 0.05, power = 0.95, alternative = "two.sided", 
    two.sided.direction = "greater", approx = FALSE, round.up = FALSE, 
    n.max = 5000, tol = 1e-07, maxiter = 1000, plot.it = TRUE, 
    add = FALSE, n.points = ifelse(x.var == "n", diff(range.x.var) + 
        1, 50), plot.col = "black", plot.lwd = 3 * par("cex"), 
    plot.lty = 1, digits = .Options$digits, ..., main = NULL, 
    xlab = NULL, ylab = NULL, type = "l") 
{
    x.var <- match.arg(x.var, c("n", "slope.over.sigma", "power", 
        "alpha"))
    y.var <- match.arg(y.var, c("power", "slope.over.sigma", 
        "n"))
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    two.sided.direction <- match.arg(two.sided.direction, c("greater", 
        "less"))
    if (x.var == y.var) 
        stop("'x.var' and 'y.var' cannot denote the same quantity")
    if (missing(range.x.var)) {
        range.x.var = switch(x.var, n = c(3, 25), slope.over.sigma = switch(alternative, 
            greater = c(0.1, 1), less = -c(1, 0.1), two.sided = if (two.sided.direction == 
                "greater") c(0.1, 1) else -c(1, 0.1)), power = c(alpha + 
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
    }, slope.over.sigma = {
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
    if (x.var != "n" && y.var != "n" && (!is.vector(n, mode = "numeric") || 
        length(n) != 1 || n != trunc(n) || n < 3)) 
        stop("'n' must be an integer greater than 2")
    if (x.var != "slope.over.sigma" && y.var != "slope.over.sigma" && 
        (!is.vector(slope.over.sigma, mode = "numeric") || length(slope.over.sigma) != 
            1)) 
        stop("'slope.over.sigma' must be a numeric scalar")
    if (x.var != "power" && y.var != "power" && (!is.vector(power, 
        mode = "numeric") || length(power) != 1 || power <= 0 || 
        power >= 1)) 
        stop("'power' must be a numeric scalar between 0 and 1")
    if (!is.vector(n.max, mode = "numeric") || length(n.max) != 
        1 || !is.finite(n.max) || n.max != trunc(n.max) || n.max < 
        3) 
        stop("'n.max' must be a positive integer greater than 2")
    if (!is.vector(maxiter, mode = "numeric") || length(maxiter) != 
        1 || !is.finite(maxiter) || maxiter != trunc(maxiter) || 
        maxiter < 2) 
        stop("'maxiter' must be a positive integer greater than 1")
    slope.string <- "slope / sigma"
    alt.string <- switch(alternative, two.sided = "(Two-Sided Alternative)", 
        less = "(Lower One-Sided Alternative)", greater = "(Upper One-Sided Alternative)")
    n.string <- paste("n =", n)
    t.string <- "t-Test for Linear Trend with\n"
    if (plot.it) 
        gen.gp.list <- checkGraphicsPars(...)$gen.gp.list
    combo <- paste(sort(c(x.var, y.var)), collapse = " & ")
    switch(combo, `alpha & slope.over.sigma` = {
        x <- seq(range.x.var[1], range.x.var[2], length = n.points)
        y <- linearTrendTestScaledMds(n = n, alpha = x, power = power, 
            alternative = alternative, two.sided.direction = two.sided.direction, 
            approx = approx, tol = tol, maxiter = maxiter)
        if (is.null(ylab)) ylab <- slope.string
        if (is.null(xlab)) xlab <- "Alpha"
        if (is.null(main)) main <- paste("Slope/Sigma vs. Alpha for ", 
            t.string, n.string, " and Power = ", format(power, 
                digits = digits), " ", alt.string, sep = "")
    }, `alpha & n` = {
        if (alternative == "greater" && slope.over.sigma <= 0) stop(paste("When alternative='greater',", 
            "'slope.over.sigma' must be positive."))
        if (alternative == "less" && slope.over.sigma >= 0) stop(paste("When alternative='less',", 
            "'slope.over.sigma' must be negative."))
        x <- seq(range.x.var[1], range.x.var[2], length = n.points)
        if (is.null(ylab)) ylab <- "Sample Size (n)"
        y <- linearTrendTestN(slope.over.sigma = slope.over.sigma, 
            alpha = x, power = power, alternative = alternative, 
            approx = approx, round.up = round.up, n.max = n.max, 
            tol = tol, maxiter = maxiter)
        if (is.null(xlab)) xlab <- "Alpha"
        if (is.null(main)) main <- paste("Sample Size vs. Alpha for ", 
            t.string, "Slope/Sigma = ", format(slope.over.sigma, 
                digits = digits), " and Power = ", format(power, 
                digits = digits), " ", alt.string, sep = "")
    }, `alpha & power` = {
        x <- seq(range.x.var[1], range.x.var[2], length = n.points)
        y <- linearTrendTestPower(n = n, slope.over.sigma = slope.over.sigma, 
            alpha = x, alternative = alternative, approx = approx)
        if (is.null(xlab)) xlab <- "Alpha"
        if (is.null(ylab)) ylab <- "Power"
        if (is.null(main)) main <- paste("Power vs. Alpha for ", 
            t.string, n.string, " and Slope/Sigma = ", format(slope.over.sigma, 
                digits = digits), " ", alt.string, sep = "")
    }, `n & slope.over.sigma` = {
        if (x.var == "slope.over.sigma") {
            if (alternative == "greater" && range.x.var[1] <= 
                0) stop(paste("When alternative='greater',", 
                "and x.var='slope.over.sigma' and y.var='n', all values of", 
                "'slope.over.sigma' must be positive,", "so the first element of 'range.x.var' must be positive."))
            if (alternative == "less" && range.x.var[2] >= 0) stop(paste("When alternative='less',", 
                "and x.var='slope.over.sigma' and y.var='n', all values of", 
                "'slope.over.sigma' must be negative,", "so the second element of 'range.x.var' must be negative."))
            x <- seq(range.x.var[1], range.x.var[2], length = n.points)
            if (is.null(ylab)) ylab <- "Sample Size (n)"
            y <- linearTrendTestN(slope.over.sigma = x, alpha = alpha, 
                power = power, alternative = alternative, approx = approx, 
                round.up = round.up, n.max = n.max, tol = tol, 
                maxiter = maxiter)
            if (is.null(xlab)) xlab <- slope.string
            if (is.null(main)) main <- paste("Sample Size vs. Slope/Sigma for ", 
                t.string, "Power = ", format(power, digits = digits), 
                " and Alpha = ", format(alpha, digits = digits), 
                " ", alt.string, sep = "")
        } else {
            x <- seq(range.x.var[1], range.x.var[2], length = n.points)
            if (is.null(xlab)) xlab <- "Sample Size (n)"
            y <- linearTrendTestScaledMds(n = x, alpha = alpha, 
                power = power, alternative = alternative, two.sided.direction = two.sided.direction, 
                approx = approx, tol = tol, maxiter = maxiter)
            if (is.null(ylab)) ylab <- slope.string
            if (is.null(main)) main <- paste("Slope/Sigma vs. Sample Size for ", 
                t.string, "Power = ", format(power, digits = digits), 
                " and Alpha = ", format(alpha, digits = digits), 
                " ", alt.string, sep = "")
        }
    }, `power & slope.over.sigma` = {
        if (x.var == "slope.over.sigma") {
            x <- seq(range.x.var[1], range.x.var[2], length = n.points)
            y <- linearTrendTestPower(n = n, slope.over.sigma = x, 
                alpha = alpha, alternative = alternative, approx = approx)
            if (is.null(xlab)) xlab <- slope.string
            if (is.null(ylab)) ylab <- "Power"
            if (is.null(main)) main <- paste("Power vs. Slope/Sigma for ", 
                t.string, n.string, " and Alpha = ", format(alpha, 
                  digits = digits), " ", alt.string, sep = "")
        } else {
            x <- seq(range.x.var[1], range.x.var[2], length = n.points)
            y <- linearTrendTestScaledMds(n = n, alpha = alpha, 
                power = x, alternative = alternative, two.sided.direction = two.sided.direction, 
                approx = approx, tol = tol, maxiter = maxiter)
            if (is.null(ylab)) ylab <- slope.string
            if (is.null(xlab)) xlab <- "Power"
            if (is.null(main)) main <- paste("Slope/Sigma vs. Power for ", 
                t.string, n.string, " and Alpha = ", format(alpha, 
                  digits = digits), " ", alt.string, sep = "")
        }
    }, `n & power` = {
        if (x.var == "n") {
            x <- seq(range.x.var[1], range.x.var[2], length = n.points)
            if (is.null(xlab)) xlab <- "Sample Size (n)"
            y <- linearTrendTestPower(n = x, slope.over.sigma = slope.over.sigma, 
                alpha = alpha, alternative = alternative, approx = approx)
            if (is.null(ylab)) ylab <- "Power"
            if (is.null(main)) main <- paste("Power vs. Sample Size for ", 
                t.string, "Slope/Sigma = ", format(slope.over.sigma, 
                  digits = digits), " and Alpha = ", format(alpha, 
                  digits = digits), " ", alt.string, sep = "")
        } else {
            if (alternative == "greater" && slope.over.sigma <= 
                0) stop(paste("When alternative='greater',", 
                "'slope.over.sigma' must be positive."))
            if (alternative == "less" && slope.over.sigma >= 
                0) stop(paste("When alternative='less',", "'slope.over.sigma' must be negative."))
            x <- seq(range.x.var[1], range.x.var[2], length = n.points)
            if (is.null(ylab)) ylab <- "Sample Size (n)"
            y <- linearTrendTestN(slope.over.sigma = slope.over.sigma, 
                alpha = alpha, power = x, alternative = alternative, 
                approx = approx, round.up = round.up, n.max = n.max, 
                tol = tol, maxiter = maxiter)
            if (is.null(xlab)) xlab <- "Power"
            if (is.null(main)) main <- paste("Sample Size vs. Power for ", 
                t.string, "Slope/Sigma = ", format(slope.over.sigma, 
                  digits = digits), " and Alpha = ", format(alpha, 
                  digits = digits), " ", alt.string)
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
