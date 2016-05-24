## A generic function for eigenvectors plots and similar graphs
vectorplot <- function (x, ...)
	UseMethod("vectorplot")

vectorplot.default <- function (x, y, col = par("col"), circle.col = "gray",
ar.length = 0.1, pos = NULL, cex = par("cex"), labels = NULL, ...)
{
	plot(x, y, type = "n", xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), asp = 1,
		...)
	abline(h = 0, col = circle.col)
	abline(v = 0, col = circle.col)
	a <- seq(0, 2 * pi, len = 100)
	lines(cos(a), sin(a), col = circle.col)
	arrows(0, 0, x, y, col = col, length = ar.length, ...)
	if (!is.null(labels)){
		## If pos is NULL, calculate pos for each variable so that label is
		## located outside
		if (is.null(pos))
			pos <- c(2, 1, 4, 3, 2)[floor((atan2(y, x)/pi + 1.25) / 0.5) + 1]
		text(x, y, labels = labels, col = col, pos = pos, cex = cex, ...)
	}
	return(invisible())	
}

vectorplot.loadings <- function (x, choices = 1L:2L, col = par("col"),
circle.col = "gray", ar.length = 0.1, pos = NULL, cex = par("cex"),
labels = rownames(x), main = deparse(substitute(x)), ...) {
	X <- x[, choices]
	vectorplot.default(X[, 1], X[, 2], col = col, circle.col = circle.col,
		ar.length = ar.length, pos = pos, cex = cex, labels = labels,
		main = main, ...)
	return(invisible(x))
}

## Plot vectors inside a circle for correlations along 2 axes (i.e., 2 columns
## in the correlation matrix). This is the typical correlations plot in PCA
vectorplot.correlation <- function (x, choices = 1L:2L, col = par("col"),
circle.col = "gray", ar.length = 0.1, pos = NULL, cex = par("cex"),
labels = rownames(x), main = deparse(substitute(x)), ...)
{
	X <- x[, choices]
	vectorplot.default(X[, 1], X[, 2], col = col, circle.col = circle.col,
		ar.length = ar.length, pos = pos, cex = cex, labels = labels,
		main = main, ...)
	return(invisible(x))
}
