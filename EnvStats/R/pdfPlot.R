pdfPlot <-
function (distribution = "norm", param.list = list(mean = 0, 
    sd = 1), left.tail.cutoff = ifelse(is.finite(supp.min), 0, 
    0.001), right.tail.cutoff = ifelse(is.finite(supp.max), 0, 
    0.001), plot.it = TRUE, add = FALSE, n.points = 1000, pdf.col = "black", 
    pdf.lwd = 3 * par("cex"), pdf.lty = 1, curve.fill = !add, 
    curve.fill.col = "cyan", x.ticks.at.all.x.max = 15, hist.col = ifelse(add, 
        "black", "cyan"), density = 5, digits = .Options$digits, 
    ..., type = "l", main = NULL, xlab = NULL, ylab = NULL, xlim = NULL, 
    ylim = NULL) 
{
    if (!is.list(param.list)) 
        stop("'param.list' must be a list.")
    check.da.list <- check.distribution.args(distribution, param.list)
    dist.abb <- check.da.list$dist.abb
    dist.name <- check.da.list$dist.name
    dist.type <- check.da.list$dist.type
    n.dist.params <- check.da.list$n.dist.params
    dist.params.names <- check.da.list$dist.params.names
    param.list <- check.da.list$param.list
    supp.min <- eval(parse(text = EnvStats::Distribution.df[dist.abb, 
        "Support.Min"]), envir = param.list)
    supp.max <- eval(parse(text = EnvStats::Distribution.df[dist.abb, 
        "Support.Max"]), envir = param.list)
    qname <- paste("q", dist.abb, sep = "")
    dname <- paste("d", dist.abb, sep = "")
    if (left.tail.cutoff == 0) {
        if (supp.min == -Inf) 
            stop(paste("The value of 'left.tail.cutoff' must be greater", 
                "than 0 for the", dist.name, "distribution since the support", 
                "on the left-hand tail is infinite."))
        else x.min <- supp.min
    }
    else x.min <- do.call(qname, c(list(p = left.tail.cutoff), 
        param.list))
    if (right.tail.cutoff == 0) {
        if (supp.max == Inf) 
            stop(paste("The value of 'right.tail.cutoff' must be greater", 
                "than 0 for the", dist.name, "distribution since the support", 
                "on the right-hand tail is infinite."))
        else x.max <- supp.max
    }
    else x.max <- do.call(qname, c(list(p = 1 - right.tail.cutoff), 
        param.list))
    if (add) {
        usr <- par("usr")
        x.min <- max(usr[1], x.min)
        x.max <- min(usr[2], x.max)
    }
    else if (!is.null(xlim)) {
        x.min <- max(xlim[1], x.min)
        x.max <- min(xlim[2], x.max)
    }
    discrete <- any(dist.type == c("Discrete", "Finite Discrete"))
    if (plot.it && !add) {
        if (is.null(xlab)) 
            xlab <- "Value of Random Variable"
        if (is.null(ylab)) 
            ylab <- ifelse(discrete, "Probability", "Relative Frequency")
        check.gp.list <- checkGraphicsPars(...)
        gp.names <- check.gp.list$gp.names
        n.gp <- check.gp.list$n.gp
        gen.gp.list <- check.gp.list$gen.gp.list
        if (is.null(main)) {
            if (any(dist.abb == c("beta", "chisq", "f", "t"))) {
                if (param.list$ncp > 0) 
                  main <- paste("Non-central ", dist.name, " Density\n", 
                    "(", paste(paste(dist.params.names, signif(unlist(param.list), 
                      digits), sep = "="), collapse = ", "), 
                    ")", sep = "")
                else {
                  main <- paste(dist.name, " Density\n", "(", 
                    paste(paste(dist.params.names[-n.dist.params], 
                      signif(unlist(param.list[-n.dist.params]), 
                        digits), sep = "="), collapse = ", "), 
                    ")", sep = "")
                }
            }
            else main <- paste(dist.name, " Density\n", "(", 
                paste(paste(dist.params.names, signif(unlist(param.list), 
                  digits), sep = "="), collapse = ", "), ")", 
                sep = "")
        }
    }
    if (dist.type == "Continuous") {
        x <- seq(x.min, x.max, len = n.points)
        y <- do.call(dname, c(list(x = x), param.list))
        if (any(dist.abb == c("beta", "exp", "gamma", "gammaAlt", 
            "weibull")) && any(index <- x == 0) && any(abs(y[index] - 
            do.call(dname, c(list(x = .Machine$double.eps), param.list))) > 
            .Machine$double.eps)) {
            x <- x[!index]
            y <- y[!index]
        }
        if (dist.abb == "beta" && any(index <- x == 1) && any(abs(y[index] - 
            do.call(dname, c(list(x = .Machine$double.eps), param.list))) > 
            .Machine$double.eps)) {
            x <- x[!index]
            y <- y[!index]
        }
        n.points <- length(x)
        if (plot.it) {
            if (!add) {
                if (is.null(xlim)) 
                  xlim <- range(x)
                if (is.null(ylim)) 
                  ylim <- c(0, max(y))
                plot(x, y, type = "n", ..., xlab = xlab, ylab = ylab, 
                  main = main, xlim = xlim, ylim = ylim)
                arg.list <- list(x = x, y = y)
                arg.list <- c(arg.list, gen.gp.list, list(type = type, 
                  col = pdf.col, lwd = pdf.lwd, lty = pdf.lty))
                do.call("lines", arg.list)
            }
            else {
                lines(x, y, ..., type = type, col = pdf.col, 
                  lwd = pdf.lwd, lty = pdf.lty)
            }
            if (curve.fill) {
                polygon(c(x, rev(x)), c(y, rep(0, n.points)), 
                  border = FALSE, col = curve.fill.col)
            }
        }
    }
    else if (discrete) {
        x <- ceiling(x.min):floor(x.max)
        y <- do.call(dname, c(list(x = x), param.list))
        nx <- length(x)
        con <- 0.4 + (0.1 * (nx - 2))/nx
        xleft <- x - con
        xright <- x + con
        ybottom <- rep(0, nx)
        if (plot.it) {
            if (!add) {
                if (is.null(xlim)) 
                  xlim.to.use <- c(min(xleft), max(xright))
                if (is.null(ylim)) 
                  ylim <- c(0, max(y))
                plot(x, y, type = "n", xaxt = "n", bty = "n", 
                  ..., xlim = xlim.to.use, ylim = ylim, xlab = xlab, 
                  ylab = ylab, main = main)
                rect(xleft = xleft, ybottom = ybottom, xright = xright, 
                  ytop = y, col = hist.col, border = pdf.col, 
                  ...)
                if (is.null(xlim) && length(x) <= x.ticks.at.all.x.max) {
                  arg.list <- list(side = 1, at = x, labels = x)
                }
                else {
                  arg.list <- list(side = 1)
                }
                arg.list <- c(arg.list, gen.gp.list)
                do.call("axis", arg.list)
            }
            else {
                o.par <- par(new = TRUE, xaxs = "d", yaxs = "d")
                on.exit(par(o.par))
                rect(xleft = xleft, ybottom = ybottom, xright = xright, 
                  ytop = y, col = hist.col, density = density, 
                  border = pdf.col, ...)
            }
        }
    }
    else {
        if (dist.name == "Zero-Modified Lognormal (Delta)" && 
            left.tail.cutoff != 0) 
            stop(paste("The value of 'left.tail.cutoff' must be 0", 
                "for the Zero-Modified Lognormal (Delta) distribution."))
        x <- seq(x.min, x.max, len = n.points)
        y <- do.call(dname, c(list(x = x), param.list))
        if (dist.name == "Zero-Modified Normal") {
            if (any(index <- x == 0)) {
                x <- c(0, x[!index])
                y <- c(y[index], y[!index])
            }
            else {
                x <- c(0, x)
                y <- c(do.call(dname, c(list(x = 0), param.list)), 
                  y)
                n.points <- n.points + 1
            }
        }
        if (plot.it) {
            if (!add) {
                if (is.null(xlim)) 
                  xlim <- range(x)
                if (is.null(ylim)) 
                  ylim <- c(0, max(y))
                plot(x[-1], y[-1], type = "n", ..., xlab = xlab, 
                  ylab = ylab, xlim = xlim, ylim = ylim, main = main)
                arg.list <- list(x = x[-1], y = y[-1])
                arg.list <- c(arg.list, gen.gp.list, list(type = type, 
                  col = pdf.col, lwd = pdf.lwd, lty = pdf.lty))
                do.call("lines", arg.list)
                if (curve.fill) {
                  polygon(c(x[-1], rev(x[-1])), c(y[-1], rep(0, 
                    n.points - 1)), border = FALSE, col = curve.fill.col)
                }
                arg.list <- list(x = x[1], y = y[1], type = "h")
                arg.list <- c(arg.list, gen.gp.list, list(col = pdf.col, 
                  lwd = pdf.lwd, lty = pdf.lty))
                do.call("points", arg.list)
            }
            else {
                lines(x[-1], y[-1], ..., type = type, col = pdf.col, 
                  lwd = pdf.lwd, lty = pdf.lty)
                if (curve.fill) {
                  polygon(c(x[-1], rev(x[-1])), c(y[-1], rep(0, 
                    n.points - 1)), border = FALSE, col = curve.fill.col)
                }
                points(x[1], y[1], type = "h", ..., col = pdf.col, 
                  lwd = pdf.lwd, lty = pdf.lty)
            }
        }
    }
    invisible(list(Quantiles = x, Probability.Densities = y))
}
