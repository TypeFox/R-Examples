plot.stepfun <- function (x, xval, xlim, ylim = range(c(y, Fn.kn)), xlab = "x", 
    ylab = "f(x)", main = NULL, add = FALSE, verticals = TRUE, 
    do.points = TRUE, pch = par("pch"), col = par("col"), col.points = col, 
    cex.points = par("cex"), col.hor = col, col.vert = col, lty = par("lty"), 
    lwd = par("lwd"), ...) 
{
    if (!is.stepfun(x)) {
        if (is.numeric(x)) {
            sarg <- substitute(x)
            x <- ecdf(x)
            attr(x, "call") <- call("ecdf", sarg)
        }
        else stop("'plot.stepfun' called with wrong type of argument 'x'")
    }
    if (missing(main)) 
        main <- {
            cl <- attr(x, "call")
            deparse(if (!is.null(cl)) 
                cl
            else sys.call())
        }
    knF <- knots(x)
    extend <- missing(xlim)
    xval <- if (missing(xval)) 
        knF
    else sort(xval)
    if (extend) {
	rx <- range(xval)
	dr <- if (length(xval) > 1) 
                max(0.08 * diff(rx), median(diff(xval)))
	else abs(xval)/16
	xlim <- rx + dr * c(-1, 1)
    }
    else {
	dr <- 0
	eps <- sqrt(.Machine$double.eps)
	xtra <- c(x(min(knF)-eps),x(max(knF)+eps))
    }
    knF <- knF[xlim[1L] - dr <= knF & knF <= xlim[2L] + dr]
    ti <- c(xlim[1L] - dr, knF, xlim[2L] + dr)
    ti.l <- ti[-length(ti)]
    ti.r <- ti[-1L]
    y <- x(0.5 * (ti.l + ti.r))
    n <- length(y)
    Fn.kn <- x(knF)
    if (add) 
        segments(ti.l, y, ti.r, y, col = col.hor, lty = lty, 
            lwd = lwd, ...)
    else {
        if (missing(ylim)) 
            ylim <- if(extend) range(y, Fn.kn) else range(y,Fn.kn,xtra)
        plot(0, 0, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, 
            ylab = ylab, main = main, ...)
        segments(ti.l, y, ti.r, y, col = col.hor, lty = lty, 
            lwd = lwd)
    }
    if (do.points) 
        points(knF, Fn.kn, pch = pch, col = col.points, cex = cex.points)
    if (verticals) 
        segments(knF, y[-n], knF, y[-1L], col = col.vert, lty = lty, 
            lwd = lwd)
    invisible(list(t = ti, y = y))
}
