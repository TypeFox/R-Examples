##
##  e r r o r b a r . R  Plot Error Bars
##


#-- Plot or add horizontal and vertical error bars
errorbar <- function(x, y, xerr = NULL, yerr = NULL,
                     bar.col = "red", bar.len = 0.01,
                     grid = TRUE, with = TRUE, add = FALSE, ...)
{
    stopifnot(is.numeric(x), is.numeric(y))
    if (!is.vector(x) || !is.vector(y) || length(x) != length(y))
        stop("Arguments 'x' and 'y' must be numeric vectors of equal length.")

    n <- length(x)
    if ( !is.null(xerr) &&
        (!is.vector(xerr) || (length(xerr) != n && length(xerr) != 1)) )
        stop("Argument 'xerr' must be a vector the same length as 'x'.")
    if ( !is.null(yerr) &&
        (!is.vector(yerr) || (length(yerr) != n && length(yerr) != 1)) )
        stop("Argument 'yerr' must be a vector the same length as 'y'.")

    if (is.null(xerr)){
        x1 <- min(x); x2 <- max(x)
    } else {
	    x1 <- min(x - abs(xerr))
	    x2 <- max(x + abs(xerr))
    }
    if (is.null(yerr)) {
        y1 <- min(y); y2 <- max(y)
    } else {
	    y1 <- min(y - abs(yerr))
	    y2 <- max(y + abs(yerr))
    }

	if (!add) {
		plot(x, y, xlim = c(x1, x2), ylim = c(y1, y2), ...)
		if (grid) grid()
	}

	# Plot the error bars
	if (!is.null(yerr))
	    segments(x, y-yerr, x, y+yerr, col = bar.col)
	if (!is.null(xerr))
	    segments(x-xerr, y, x+xerr, y, col = bar.col)


	if (with) {
		xd <- bar.len * (x2 - x1) / 2.0
		yd <- bar.len * (y2 - y1) / 2.0

        if (!is.null(yerr)) {
		    segments(x-xd, y-yerr, x+xd, y-yerr, col = bar.col)
		    segments(x-xd, y+yerr, x+xd, y+yerr, col = bar.col)
	    }
	    if (!is.null(xerr)) {
		    segments(x-xerr, y-yd, x-xerr, y+yd, col = bar.col)
		    segments(x+xerr, y-yd, x+xerr, y+yd, col = bar.col)
	    }
	}

	invisible()
}
