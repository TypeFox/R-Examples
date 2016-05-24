errorBar <-
function (x, y = NULL, lower, upper, incr = TRUE, draw.lower = TRUE, 
    draw.upper = TRUE, bar.ends = TRUE, gap = TRUE, add = FALSE, 
    horizontal = FALSE, gap.size = 0.75, bar.ends.size = 1, col = 1, 
    ..., xlab = deparse(substitute(x)), xlim, ylim) 
{
    draw.null.warn <- function(draw, gap) {
        if (!any(draw)) {
            warning("Not enough room for a gap.")
            draw <- !draw
            gap <- 0
        }
        invisible(list(draw = draw, gap = gap))
    }
    if (missing(x)) 
        stop("no data for x or y")
    if (missing(y)) {
        if (missing(xlab)) 
            xlab <- "Index"
        y <- x
        x <- time(x)
    }
    n <- length(x)
    if (length(y) != n) 
        stop("length of y must equal the length of x")
    center <- if (horizontal) 
        x
    else y
    if (missing(lower)) 
        stop("you must provide lower")
    if (length(lower) > 1 && length(lower) != n) 
        stop("length of lower must be 1 or equal to the length of x")
    if (incr) 
        lower <- center - abs(lower)
    else lower <- rep(lower, length = n)
    if (draw.lower && any(lower >= center)) 
        warning(paste("There are values of 'lower' which are greater or equal to ", 
            if (horizontal) 
                "x"
            else "y"))
    if (missing(upper)) 
        upper <- 2 * center - lower
    else {
        if (length(upper) > 1 && length(upper) != n) 
            stop("length of upper must be 1 or equal to the length of x")
        if (incr) 
            upper <- center + upper
        else upper <- rep(upper, length = n)
    }
    if (draw.upper && any(upper <= center)) 
        warning(paste("There are values of 'upper' which are smaller or\nequal to ", 
            if (horizontal) 
                "x"
            else "y"))
    if (!add) {
        if (horizontal) {
            if (missing(ylim)) 
                plot(x, y, xlim = if (missing(xlim)) 
                  range(c(lower, upper), na.rm = TRUE)
                else xlim, xlab = xlab, col = col, ...)
            else plot(x, y, xlim = if (missing(xlim)) 
                range(c(lower, upper), na.rm = TRUE)
            else xlim, ylim = ylim, xlab = xlab, col = col, ...)
        }
        else {
            if (missing(xlim)) 
                plot(x, y, ylim = if (missing(ylim)) 
                  range(c(lower, upper), na.rm = TRUE)
                else ylim, xlab = xlab, col = col, ...)
            else plot(x, y, ylim = if (missing(ylim)) 
                range(c(lower, upper), na.rm = TRUE)
            else ylim, xlim = xlim, xlab = xlab, col = col, ...)
        }
    }
    col <- rep(col, length = n)
    check.gp.list <- checkGraphicsPars(...)
    gen.gp.list <- check.gp.list$gen.gp.list
    if (horizontal) {
        if (gap) 
            gap <- gap.size * par("cxy")[1]
        if (bar.ends) 
            size.bar <- bar.ends.size * par("cxy")[2]
        if (draw.lower) {
            draw <- x - lower > gap
            z <- draw.null.warn(draw, gap)
            draw <- z$draw
            gap <- z$gap
            arg.list <- c(list(x0 = lower[draw], y0 = y[draw], 
                x1 = x[draw] - gap, y1 = y[draw], col = col[draw]), 
                gen.gp.list)
            do.call(segments, arg.list)
            if (bar.ends) {
                arg.list <- c(list(x0 = lower[draw], y0 = y[draw] - 
                  size.bar, x1 = lower[draw], y1 = y[draw] + 
                  size.bar, col = col[draw]), gen.gp.list)
                do.call(segments, arg.list)
            }
        }
        if (draw.upper) {
            draw <- upper - x > gap
            z <- draw.null.warn(draw, gap)
            draw <- z$draw
            gap <- z$gap
            arg.list <- c(list(x0 = x[draw] + gap, y0 = y[draw], 
                x1 = upper[draw], y1 = y[draw], col = col[draw]), 
                gen.gp.list)
            do.call(segments, arg.list)
            if (bar.ends) {
                arg.list <- c(list(x0 = upper[draw], y0 = y[draw] - 
                  size.bar, x1 = upper[draw], y1 = y[draw] + 
                  size.bar, col = col[draw]), gen.gp.list)
                do.call(segments, arg.list)
            }
        }
        ret.list <- list(group.centers = y, group.stats = cbind(Center = x, 
            Lower = lower, Upper = upper))
    }
    else {
        if (gap) 
            gap <- gap.size * par("cxy")[2]
        if (bar.ends) 
            size.bar <- bar.ends.size * par("cxy")[1]
        if (draw.upper) {
            draw <- upper - y > gap
            z <- draw.null.warn(draw, gap)
            draw <- z$draw
            gap <- z$gap
            arg.list <- c(list(x0 = x[draw], y0 = y[draw] + gap, 
                x1 = x[draw], y1 = upper[draw], col = col[draw]), 
                gen.gp.list)
            do.call(segments, arg.list)
            if (bar.ends) {
                arg.list <- c(list(x0 = x[draw] - size.bar, y0 = upper[draw], 
                  x1 = x[draw] + size.bar, y1 = upper[draw], 
                  col = col[draw]), gen.gp.list)
                do.call(segments, arg.list)
            }
        }
        if (draw.lower) {
            draw <- y - lower > gap
            z <- draw.null.warn(draw, gap)
            draw <- z$draw
            gap <- z$gap
            arg.list <- c(list(x0 = x[draw], y0 = y[draw] - gap, 
                x1 = x[draw], y1 = lower[draw], col = col[draw]), 
                gen.gp.list)
            do.call(segments, arg.list)
            if (bar.ends) {
                arg.list <- c(list(x0 = x[draw] - size.bar, y0 = lower[draw], 
                  x1 = x[draw] + size.bar, y1 = lower[draw], 
                  col = col[draw]), gen.gp.list)
                do.call(segments, arg.list)
            }
        }
        ret.list <- list(group.centers = x, group.stats = cbind(Center = y, 
            Lower = lower, Upper = upper))
    }
    invisible(ret.list)
}
