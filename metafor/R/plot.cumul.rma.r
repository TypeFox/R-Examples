plot.cumul.rma <-
function (x, yaxis = "tau2", xlim, ylim, xlab, ylab, at, transf, 
    atransf, targs, digits, cols = c("gray90", "gray10"), addgrid = TRUE, 
    pch = 19, cex = 1, lwd = 2, ...) 
{
    if (!is.element("cumul.rma", class(x))) 
        stop("Argument 'x' must be an object of class \"cumul.rma\".")
    if (is.null(x$tau2)) 
        stop("Either a fixed-effects model or not sufficient data to estimate tau^2.")
    yaxis <- match.arg(yaxis, c("tau2", "I2", "H2"))
    if (missing(transf)) 
        transf <- FALSE
    if (missing(atransf)) 
        atransf <- FALSE
    transf.char <- deparse(substitute(transf))
    atransf.char <- deparse(substitute(atransf))
    if (transf.char != "FALSE" && atransf.char != "FALSE") 
        stop("Use either 'transf' or 'atransf' to specify a transformation (not both).")
    if (missing(xlab)) 
        xlab <- .setlab(x$measure, transf.char, atransf.char, 
            gentype = 2)
    if (missing(ylab)) {
        if (yaxis == "tau2") 
            ylab <- "Amount of Heterogeneity (tau^2)"
        if (yaxis == "I2") 
            ylab <- "Percentage of Variability due to Heterogeneity (I^2)"
        if (yaxis == "H2") 
            ylab <- "Ratio of Total Variability to Sampling Variability (H^2)"
    }
    if (missing(at)) 
        at <- NULL
    if (missing(targs)) 
        targs <- NULL
    if (missing(digits)) {
        if (yaxis == "tau2") 
            digits <- c(2, 3)
        if (yaxis == "I2") 
            digits <- c(2, 1)
        if (yaxis == "H2") 
            digits <- c(2, 1)
    }
    else {
        if (length(digits) == 1L) 
            digits <- c(digits, digits)
    }
    dat <- data.frame(estim = x$estimate)
    if (yaxis == "tau2") 
        dat$yval <- x$tau2
    if (yaxis == "I2") 
        dat$yval <- x$I2
    if (yaxis == "H2") 
        dat$yval <- x$H2
    dat <- na.omit(dat)
    k <- nrow(dat)
    if (is.function(transf)) {
        if (is.null(targs)) {
            dat$estim <- sapply(dat$estim, transf)
        }
        else {
            dat$estim <- sapply(dat$estim, transf, targs)
        }
    }
    if (missing(xlim)) {
        xlim <- range(dat$estim)
    }
    else {
        xlim <- sort(xlim)
    }
    if (missing(ylim)) {
        ylim <- range(dat$yval)
    }
    else {
        ylim <- sort(ylim)
    }
    if (!is.null(at)) {
        xlim[1] <- min(c(xlim[1], at), na.rm = TRUE)
        xlim[2] <- max(c(xlim[2], at), na.rm = TRUE)
    }
    plot(NA, NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, 
        xaxt = "n", yaxt = "n")
    if (is.null(at)) {
        at <- axTicks(side = 1)
    }
    else {
        at <- at[at > par("usr")[1]]
        at <- at[at < par("usr")[2]]
    }
    at.lab <- at
    if (is.function(atransf)) {
        if (is.null(targs)) {
            at.lab <- formatC(sapply(at.lab, atransf), digits = digits[1], 
                format = "f")
        }
        else {
            at.lab <- formatC(sapply(at.lab, atransf, targs), 
                digits = digits[1], format = "f")
        }
    }
    else {
        at.lab <- formatC(at.lab, digits = digits[1], format = "f")
    }
    axis(side = 1, at = at, labels = at.lab)
    aty <- axTicks(side = 2)
    axis(side = 2, at = aty, labels = formatC(aty, digits = digits[2], 
        format = "f"))
    if (addgrid) {
        abline(v = at, lty = "dotted", col = "lightgray")
        abline(h = aty, lty = "dotted", col = "lightgray")
    }
    cols.points <- colorRampPalette(cols)(k)
    for (i in 1:(k - 1)) {
        estims <- approx(c(dat$estim[i], dat$estim[i + 1]), n = 50)$y
        yvals <- approx(c(dat$yval[i], dat$yval[i + 1]), n = 50)$y
        cols.lines <- colorRampPalette(c(cols.points[i], cols.points[i + 
            1]))(50)
        segments(estims[-50], yvals[-50], estims[-1], yvals[-1], 
            col = cols.lines, lwd = lwd)
    }
    points(dat$estim, dat$yval, pch = pch, col = cols.points, 
        cex = cex)
    invisible(dat)
}
