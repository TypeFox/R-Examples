## A wrapper around cor() and the like, building a "correlation" S3 object
## TODO: cov.wt(), cov2correlation(), and perhaps, functions cov.XXX() from MASS
## TODO: max, min, range, which.max, which.min for 'correlation' objects that do
## not consider elements on the diagonal... or put something else to avoid it is
## extracted for max, or which.max??? + something like 'highest' which considers
## the absolute value??? How to deal with that?

## A generic function to calculate correlation from an object
correlation <- function (x, ...)
	UseMethod("correlation")

correlation.formula <- function (formula, data = NULL, subset, na.action, ...)
{
	mt <- terms(formula, data = data)
    if (attr(mt, "response") > 0L) 
        stop("response not allowed in formula")
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    mf$... <- NULL
    mf[[1L]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    
	## Copied from stats:::.check_vars_numeric()
	.check_vars_numeric <- function (mf) {
	    mt <- attr(mf, "terms")
	    mterms <- attr(mt, "factors")
	    mterms <- rownames(mterms)[apply(mterms, 1L, function(x) any(x > 
	        0L))]
	    any(sapply(mterms, function(x) is.factor(mf[, x]) || !is.numeric(mf[, 
	        x])))
	}

	if (.check_vars_numeric(mf)) 
        stop("Correlation applies only to numerical variables")
    mt <- attr(mf, "terms")
    attr(mt, "intercept") <- 0L
    x <- model.matrix(mt, mf)
    res <- correlation.default(x, ...)
    cl[[1L]] <- as.name("correlation")
    attr(res, "call") <- cl
    attr(res, "na.method") <- NULL
	if (!is.null(na.action))
        attr(res, "na.action") <- as.character(na.action)
    return(res)
}

## Create the 'correlation' object (same arguments as cor() in stats package)
correlation.default <- function (x, y = NULL, use = "everything",
method = c("pearson", "kendall", "spearman"), ...)
{
	Call <- match.call()
	x <- as.matrix(x)
	na.methods <- c("all.obs", "complete.obs", "pairwise.complete.obs",
		"everything", "na.or.complete")
	na.method <- pmatch(use, na.methods)
	method <- match.arg(method)
	
	## Just call cor in stats package
	res <- stats::cor(x = x, y = y, use = use, method = method)
	
	## We want to return a correlation matrix, even if there is one correlation
	if (length(res) == 1) {
		## Create a simple correlation matrix using 'x' and 'y' as labels
		res <- matrix(c(1, res, res, 1), ncol = 2,
			dimnames = list(c("x", "y"), c("x", "y")))
	}
	
	## Same strings as for cor.test()
	attr(res, "method") <- switch(method,
		pearson = "Pearson's product-moment correlation",
		kendall = "Kendall's rank correlation tau",
		spearman = "Spearman's rank correlation rho",
		method)
	attr(res, "na.method") <- na.methods[na.method]
	attr(res, "call") <- Call
	class(res) <- c("correlation", "matrix")
	
	return(res)
}

## Check if an object is a correlation matrix
is.correlation <- function (x)
	return(inherits(x, "correlation"))

## Transform a square matrix or a data.frame with values between -1 and 1
## in a 'correlation' object
## TODO: should we keep more attributes, in order to document other correlation
## calculations?
as.correlation <- function (x) {
	if (is.correlation(x)) return(x)
	
	## Make sure we have a matrix with numeric data, dimnames and nothing else
	## (drop all other arguments, except 'comment', perhaps)
	res <- structure(as.numeric(x), dim = dim(x), dimnames = dimnames(x))
	
	## Check that it is a square (2D) matrix, or an atomic number
	d <- dim(x)
	if (is.null(d)) {
		## Is this an atomic number?
		if (length(x) == 1) {
			## Create the simplest correlation matrix using
			## generic 'x' and 'y' labels
			res <- matrix(c(1, res, res, 1), ncol = 2,
				dimnames = list(c("x", "y"), c("x", "y")))
		}
	} else {  # Check that it is a square matrix
	if (length(d) != 2 || d[1] != d[2])
		stop("x must be a square matrix")
	}
	
	## Check the range that must be between -1 and 1
	rg <- range(res, na.rm = TRUE)
	if (rg[1] < -1 || rg[2] > 1)
		stop("A correlation matrix cannot have values lower than -1 or larger than 1")
	
	## Reinject comment, if it exists
	comment(res) <- comment(x)
	
	## Look for a "method" attribute
	attr(res, "method") <- attr(x, "method")
	## ... and a na.method, or na.action attribute
	attr(res, "na.action") <- attr(x, "na.action")
	attr(res, "na.method") <- attr(x, "na.method")
	
	## Set this as both a 'correlation' and 'matrix' S3 object
	class(res) <- c("correlation", "matrix")
	
	return(res)
}

## Print a 'correlation' object
print.correlation <- function (x, digits = 3, cutoff = 0, ...)
{
	if (!is.correlation(x))
		stop("x must be a 'correlation' object (correlation matrix)")
	
	method <- attr(x, "method")
	if (is.null(method)) {
		cat("Correlation matrix:\n")
	} else {
		cat("Matrix of ", method, ":\n", sep = "")
	}
	
	na.method <- attr(x, "na.method")
	if (!is.null(na.method)) {
		cat("(calculation uses ", na.method, ")\n", sep = "")
	} else {
		na.action <- attr(x, "na.action")
		if (!is.null(na.action))
			cat("(missing values are managed with ", na.action, ")\n", sep = "")
	}
	cform <- format(round(x, digits = digits))
	nc <- nchar(cform[1L], type = "c")
	cform[abs(x) < cutoff] <- paste(rep(" ", nc), collapse = "")
	print(cform, quote = FALSE, ...)
	
	return(invisible(x))	
}

## Summary of a 'correlation' object
summary.correlation <- function (object, cutpoints = c(0.3, 0.6, 0.8, 0.9, 0.95),
symbols = c(" ", ".", ",", "+", "*", "B"), ...)
{
	## Replace the correlation matrix by symbols using symnum()
	res <- symnum(unclass(object), cutpoints = cutpoints, symbols = symbols,
		corr = TRUE, ...)
	
	## Reinject comment, if it exists
	comment(res) <- comment(object)
	
	## Look for a "method" attribute
	attr(res, "method") <- attr(object, "method")
	## ... and na.action/na.method attributes
	attr(res, "na.action") <- attr(object, "na.action")
	attr(res, "na.method") <- attr(object, "na.method")
	
	## Set this as 'summary.correlation' object
	class(res) <- c("summary.correlation", "noquote")
	
	return(res)
}

## Print method for the 'summary.correlation' object
print.summary.correlation <- function (x, ...)
{	
	method <- attr(x, "method")
	if (is.null(method)) {
		cat("Correlation matrix:\n")
	} else {
		cat("Matrix of ", method, ":\n", sep = "")
	}
	
	na.method <- attr(x, "na.method")
	if (!is.null(na.method)) {
		cat("(calculation uses ", na.method, ")\n", sep = "")
	} else {
		na.action <- attr(x, "na.action")
		if (!is.null(na.action))
			cat("(missing values are managed with ", na.action, ")\n", sep = "")
	}
	
	print(structure(as.character(x), dim = dim(x), dimnames = dimnames(x),
		legend = attr(x, "legend"), class = "noquote"), ...)
	
	return(invisible(x))	
}

## Plot a 'correlation' object (basically the ellipse's plotcorr() function, but
## as plot() method for 'corr' object and with different default values
## Also, numbers are printed inside the ellipses with numbers = TRUE
## TODO: change the way labels are plotted
## TODO: a comparison plot, when y is not NULL
plot.correlation <- function (x, y = NULL, outline = TRUE,
cutpoints = c(0.3, 0.6, 0.8, 0.9, 0.95), palette = rwb.colors, col = NULL,
numbers = TRUE, digits = 2, type = c("full", "lower", "upper"),
diag = (type == "full"), cex.lab = par("cex.lab"), cex = 0.75 * par("cex"), ...)
{
	if (!is.correlation(x))
		stop("x must be a 'correlation' object")
	
	type <- match.arg(type)
	diag <- as.logical(diag[1])
	## Compute colors from cutpoints and palette
	if (is.null(col)) {
		## -1.1 to include -1 - intervals are (,]
		## cutpoints - 0.0001 for positive values to include lower limit instead
		br <- c(-1.1, rev(-cutpoints), cutpoints - 0.0001, 1)
		ct <- cut(x, breaks = br)
		col <- palette(length(levels(ct)))[as.numeric(ct)]
	}
	
	## Call the plotcorr() function from ellipse package
	plotcorr(x, outline = outline, col = col, numbers = FALSE, type = type,
		diag = diag, cex.lab = cex.lab, cex = cex, ...)
	## Do we print the numbers inside the ellipses?
	if (isTRUE(numbers)) {
		coords <- expand.grid(1:nrow(x), nrow(x):1)
		labels <- format(round(x, digits = digits), digits = digits)
		## Do we plotted only upper or lower triangle and diagonal?
		## Note: we need to invert y-coordinates!
		yinv <- max(coords) + 1 - coords[, 2] 
		if (diag) {
			if (type == "lower") {
				## Keep only lower triangle + diagonal
				coords <- coords[coords[, 1] <= yinv, ]
				coords <- coords[order(coords[, 1]), ]
				labels <- labels[lower.tri(labels, diag = TRUE)]
			} else if (type == "upper") {
				## Keep only upper triangle
				coords <- coords[coords[, 1] >= yinv, ]
				coords <- coords[order(coords[, 1]), ]
				labels <- labels[upper.tri(labels, diag = TRUE)]
			}
		} else {  # No diagonals
			if (type == "lower") {
				## Keep only lower triangle
				coords <- coords[coords[, 1] < yinv, ]
				coords <- coords[order(coords[, 1]), ]
				labels <- labels[lower.tri(labels)]
			} else if (type == "upper") {
				## Keep only upper triangle
				coords <- coords[coords[, 1] > yinv - 1, ]
				coords <- coords[order(coords[, 1]), ]
				coords[, 2] <- coords[, 2] - 1
				labels <- labels[upper.tri(labels)]
			} else {
				## Plot everything, except diagonal => put test to "" there
				diag(labels) <- ""
			}
		}
		text(coords, labels = labels, cex = cex, ...)
	}
	return(invisible())
}

## Add vectors for supplementary variables in a PCA correlation plot
lines.correlation <- function (x, choices = 1L:2L, col = par("col"), lty = 2,
ar.length = 0.1, pos = NULL, cex = par("cex"), labels = rownames(x),  ...)
{
	corrs <- x[, choices]
	arrows(0, 0, corrs[, 1], corrs[, 2], col = col, lty = lty,
		length = ar.length, ...)
	if (!is.null(labels)){
		## If pos is NULL, calculate pos for each variable so that label is
		## located outside
		if (is.null(pos))
			pos <- c(2, 1, 4, 3, 2)[floor((atan2(corrs[, 2], corrs[, 1])/pi +
				1.25) / 0.5) + 1]
		text(corrs, labels = labels, col = col, pos = pos, cex = cex, ...)
	}
	return(invisible(x))
}
