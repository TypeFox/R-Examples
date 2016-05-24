plotAovDesign <-
function (x.var = "n", y.var = "power", range.x.var = NULL, n.vec = c(25, 
    25), mu.vec = c(0, 1), sigma = 1, alpha = 0.05, power = 0.95, 
    round.up = FALSE, n.max = 5000, tol = 1e-07, maxiter = 1000, 
    plot.it = TRUE, add = FALSE, n.points = 50, plot.col = 1, 
    plot.lwd = 3 * par("cex"), plot.lty = 1, digits = .Options$digits, 
    main = NULL, xlab = NULL, ylab = NULL, type = "l", ...) 
{
    x.var <- match.arg(x.var, c("n", "power", "alpha"))
    y.var <- match.arg(y.var, c("power", "n"))
    if (x.var == y.var) 
        stop("'x.var' and 'y.var' cannot denote the same quantity")
    if (missing(range.x.var)) {
        range.x.var <- switch(x.var, n = c(2, 50), power = c(alpha + 
            .Machine$double.eps, 0.95), alpha = c(0.01, 0.2))
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
    if (!is.vector(mu.vec, mode = "numeric")) 
        stop("'mu.vec' must be a numeric vector")
    if (!all(is.finite(mu.vec))) 
        stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
            "Undefined (Nan) values are not allowed in", "'mu.vec'"))
    if (length(mu.vec) < 2 || length(unique(mu.vec)) == 1) 
        stop("'mu.vec' must have at least 2 distinct values")
    n.grps <- length(mu.vec)
    if (!is.vector(sigma, mode = "numeric") || length(sigma) != 
        1 || !is.finite(sigma) || sigma <= 0) 
        stop("'sigma' must be a positive number")
    if (x.var != "n" && y.var != "n") {
        if (length(n.vec) != n.grps) 
            stop("'n.vec' must be the same length as 'mu.vec'")
        if (!is.vector(n.vec, mode = "numeric")) 
            stop("'n.vec' must be a numeric vector.")
        if (!all(is.finite(n.vec))) 
            stop(paste("Missing (NA), Infinite (Inf, -Inf), and", 
                "Undefined (Nan) values are not allowed in", 
                "'n.vec'"))
        if (!all(n.vec == trunc(n.vec)) || any(n.vec < 2)) 
            stop(paste("All values of 'n.vec' must be integers", 
                "greater than or equal to 2"))
    }
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
    if (x.var != "n" && y.var != "n") {
        n.string <- paste("n = (", paste(n.vec, collapse = ","), 
            ")", sep = "")
    }
    if (x.var != "power" && y.var != "power") {
        power.string <- paste("power =", format(power, digits = digits))
    }
    if (x.var != "alpha") {
        alpha.string <- paste("alpha =", format(alpha, digits = digits))
    }
    mu.string <- paste(format(mu.vec, digits = digits), collapse = ",")
    index <- nchar(mu.string)
    last.blank <- substring(mu.string, index, index) == " "
    while (last.blank) {
        index <- index - 1
        mu.string <- substring(mu.string, 1, index)
        last.blank <- substring(mu.string, index, index) == " "
    }
    mu.string <- paste("mu = (", mu.string, ")", sep = "")
    sigma.string <- paste("sigma =", format(sigma, digits = digits))
    if (plot.it) 
        gen.gp.list <- checkGraphicsPars(...)$gen.gp.list
    combo <- paste(c(x.var, y.var), collapse = " & ")
    switch(combo, `alpha & n` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(ylab)) ylab <- "Sample Size per Group (n)"
        y <- aovN(mu.vec = mu.vec, sigma = sigma, alpha = x, 
            power = power, round.up = round.up, n.max = n.max, 
            tol = tol, maxiter = maxiter)
        if (is.null(xlab)) xlab <- "Alpha"
        if (is.null(main)) main <- paste("Sample Size vs. Alpha for One-Way ANOVA\n", 
            "with ", mu.string, ", ", sigma.string, ", and ", 
            power.string, sep = "")
    }, `alpha & power` = {
        x <- seq(min.x, max.x, length = n.points)
        y <- aovPower(n.vec = n.vec, mu.vec = mu.vec, sigma = sigma, 
            alpha = x)
        if (is.null(xlab)) xlab <- "Alpha"
        if (is.null(ylab)) ylab <- "Power"
        if (is.null(main)) main <- paste("Power vs. Alpha for One-Way ANOVA\n", 
            "with ", mu.string, ", ", n.string, ", and ", sigma.string, 
            sep = "")
    }, `n & power` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(xlab)) xlab <- "Group Sample Size (n)"
        N <- length(x)
        y <- numeric(N)
        for (i in 1:N) y[i] <- aovPower(n.vec = rep(x[i], n.grps), 
            mu.vec = mu.vec, sigma = sigma, alpha = alpha)
        if (is.null(ylab)) ylab <- "Power"
        if (is.null(main)) main <- paste("Power vs. Sample Size for One-Way ANOVA\n", 
            "with ", mu.string, ", ", sigma.string, ", and ", 
            alpha.string, sep = "")
    }, `power & n` = {
        x <- seq(min.x, max.x, length = n.points)
        if (is.null(xlab)) xlab <- "Power"
        y <- aovN(mu.vec = mu.vec, sigma = sigma, alpha = alpha, 
            power = x, round.up = round.up, n.max = n.max, tol = tol, 
            maxiter = maxiter)
        if (is.null(ylab)) ylab <- "Group Sample Size (n)"
        if (is.null(main)) main <- paste("Sample Size vs. Power for One-Way ANOVA\n", 
            "with ", mu.string, ", ", sigma.string, ", and ", 
            alpha.string, sep = "")
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
