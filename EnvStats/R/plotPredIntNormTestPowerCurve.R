plotPredIntNormTestPowerCurve <-
function (n = 8, df = n - 1, n.mean = 1, k = 1, range.delta.over.sigma = c(0, 
    5), pi.type = "upper", conf.level = 0.95, plot.it = TRUE, 
    add = FALSE, n.points = 20, plot.col = "black", plot.lwd = 3 * 
        par("cex"), plot.lty = 1, digits = .Options$digits, ..., 
    main = NULL, xlab = NULL, ylab = NULL, type = "l") 
{
    pi.type <- match.arg(pi.type, c("upper", "lower"))
    if (is.null(range.delta.over.sigma) || !all(is.finite(range.delta.over.sigma)) || 
        !is.vector(range.delta.over.sigma, mode = "numeric") || 
        length(range.delta.over.sigma) != 2) 
        stop(paste("'range.delta.over.sigma' must be a numeric vector of length 2", 
            "with no missing (NA), infinite(-Inf, Inf), or undefined(NaN) values"))
    min.x <- range.delta.over.sigma[1]
    max.x <- range.delta.over.sigma[2]
    if (min.x >= max.x) 
        stop("The second element of 'range.delta.over.sigma' must be larger than the first")
    if (!is.vector(n.points, mode = "numeric") || length(n.points) != 
        1 || n.points != trunc(n.points) || n.points < 2) 
        stop("'n.points' must be an integer larger than 1")
    if (is.null(n) || !is.finite(n) || !is.vector(n, mode = "numeric") || 
        length(n) != 1 || n < 2 || n != trunc(n)) 
        stop("'n' must be an integer greater than 1")
    if (is.null(df) || !is.finite(df) || !is.vector(df, mode = "numeric") || 
        length(df) != 1 || df != trunc(df) || any(df < 1) || 
        is.null(k) || !is.finite(k) || !is.vector(k, mode = "numeric") || 
        length(k) != 1 || k != trunc(k) || any(k < 1) || is.null(n.mean) || 
        !is.finite(n.mean) || !is.vector(n.mean, mode = "numeric") || 
        length(n.mean) != 1 || n.mean != trunc(n.mean) || any(n.mean < 
        1)) 
        stop("'df', 'k' and 'n.mean' must be positive integers")
    if (is.null(conf.level) || !is.finite(conf.level) || !is.vector(conf.level, 
        mode = "numeric") || length(conf.level) != 1 || conf.level <= 
        .Machine$double.eps || conf.level >= 1 - .Machine$double.eps) {
        stop("'conf.level' must be a scalar between 0 and 1")
    }
    delta.string <- "(Future Mean - Background Mean) / SD"
    pi.string <- switch(pi.type, lower = "(Lower One-Sided PI)", 
        upper = "(Upper One-Sided PI)")
    n.string <- paste("n =", n)
    if (df != n - 1) 
        n.string <- paste(n.string, ", df = ", df, sep = "")
    if (n.mean != 1) 
        n.string <- paste(n.string, ", n.mean = ", n.mean, sep = "")
    k.string <- paste("k =", k)
    if (plot.it) 
        gen.gp.list <- checkGraphicsPars(...)$gen.gp.list
    x <- seq(min.x, max.x, length = n.points)
    y <- predIntNormTestPower(n = n, df = df, n.mean = n.mean, 
        k = k, delta.over.sigma = x, pi.type = pi.type, conf.level = conf.level)
    if (is.null(xlab)) 
        xlab <- delta.string
    if (is.null(ylab)) 
        ylab <- "Power"
    if (is.null(main)) 
        main <- paste("Power vs. Delta/Sigma for Normal Prediction Interval with\n", 
            n.string, ", ", k.string, ", and Confidence Level = ", 
            format(conf.level, digits = digits), " ", pi.string, 
            sep = "")
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
    names(ret.list) <- c("delta.over.sigma", "power")
    invisible(ret.list)
}
