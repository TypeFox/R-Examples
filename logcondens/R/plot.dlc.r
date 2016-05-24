plot.dlc <- function(x, which = c("density", "log-density", "CDF"), add.title = TRUE, legend.pos = "topright", ...){

    object <- x
    which <- match.arg(which)
    x <- object$x
    n <- length(x)
    xn <- object$xn
    knot <- x[object$IsKnot == 1]
    if (legend.pos[1] == "none"){legend.pos <- c(Inf, Inf)}

    # compute function values
    d <- diff(range(x))
    xs <- seq(x[1], x[n], length.out = 500)
    f.mat <- matrix(NA, nrow = length(xs), ncol = 5)
    for (i in 1:length(xs)){f.mat[i, 1:3] <- evaluateLogConDens(xs[i], object, which = 1:3)[, c("log-density", "density", "CDF")]}

    xli <- range(xs) + 0.1 * c(-1, 1) * d
    yli <- range(c(0, f.mat[, 2]))
    yli.log <- range(f.mat[, 1])
    cho <- 1
    d <- diff(range(x))
    if (object$smoothed == TRUE){
        xli <- range(c(xli, object$xs))
        yli <- range(c(yli, object$f.smoothed))
        yli.log <- c(min(c(f.mat[, 1], log(object$f.smoothed))), log(max(yli)))
        cho <- 1:2
    }

    ## organize title
    if (add.title == FALSE){main.tit <- ""; o <- 1}
    if (add.title == TRUE){main.tit <- "Log-concave density estimation from i.i.d. data"; o <- 4}

    par(las = 1, mar = c(4.5, 4, o, 1))

    if (which == "density"){
        plot(0, 0, type = 'n', xlim = xli, ylim = yli, xlab = "x", ylab = "density", main = main.tit)
        rug(x)
        lines(xs, f.mat[, 2], col = 2, lwd = 2)
        if (object$smoothed == TRUE){lines(object$xs, object$f.smoothed, col = 3, lwd = 2)}
        if (add.title == TRUE){mtext("Density estimate. Dashed vertical lines indicate knots of the log-density.", 3, 0.5, cex = .75)}
        legend(legend.pos, c(expression("log-concave "*hat(f)[n]), expression("log-concave smoothed "*hat(f)[n]*"*"))[cho],
            lty = 1, lwd = 2, col = c(2, 3)[cho], bty = "n")
        segments(c(min(xli) - 10 * d, x[n]), c(0, 0), c(x[1], max(xli) + 10 * d), c(0, 0), col = 2, lwd = 2)
        segments(knot, 0, knot, 0.2 * diff(yli), lty = 2)
        }

    if (which == "log-density"){
        plot(0, 0, type = 'n', xlim = xli, ylim = yli.log, xlab = "x", ylab = "log-density", main = main.tit)
        rug(x)
        lines(xs, f.mat[, 1], col = 2, lwd = 2)
        if (object$smoothed == TRUE){lines(object$xs, log(object$f.smoothed), col = 3, lwd = 2)}
        if (add.title == TRUE){mtext("Estimate of log-density. Dashed vertical lines indicate knots of the log-density.", 3, 0.5, cex = .75)}
        legend(legend.pos, c(expression("log-concave "*hat(phi)[n]), expression("log-concave smoothed "*hat(phi)[n]*"*"))[cho],
            lty = 1, lwd = 2, col = c(2, 3)[cho], bty = "n")
        segments(knot, min(yli.log), knot, min(yli.log) + 0.2 * diff(yli.log), lty = 2)
        }

    if (which == "CDF"){
        n0 <- length(xn)
        ED <- 0 : (n0 + 1) / n0
        plot(xs, f.mat[, 3], type = "n", ylab = "Distribution functions", main = main.tit, xlim = xli)
        rug(x); lines(xs, f.mat[, 3], col = 2, lwd = 2)
        lines(c(min(xn) - 10, xn, max(xn) + 10), ED, type = 's', lwd = 1.5)
        abline(v = knot, lty = 3)
        if (object$smoothed == TRUE){lines(object$xs, object$F.smoothed, col = 3, lwd = 2)}    
        if (add.title == TRUE){mtext("Estimate of CDF based on density.", 3, 0.5, cex = .75)}
        legend(legend.pos, c("ECDF", expression("log-concave "*hat(F)[n]), expression("log-concave smoothed "*hat(F)[n]*"*"))[c(1, cho + 1)],
            lty = 1, lwd = 2, col = c(1, 2, 3)[c(1, cho + 1)], bty = "n")
        segments(c(min(xli) - 10 * d, x[n]), c(0, 1), c(x[1], max(xli) + 10 * d), c(0, 1), col = 2, lwd = 2)
        }
}
