## Define a "pcomp" S3 object for PCA, because there is too much chaos with
## default "prcomp" and "princomp" R objects, plus "pca" in ade4 and labdsv,
## "PCA" in FactoMineR, etc.

## Create the pcomp generic function that returns a "pcomp" object
pcomp <- function (x, ...)
	UseMethod("pcomp")
	
pcomp.formula <- function (formula, data = NULL, subset, na.action,
method = c("svd", "eigen"), ...)
{
	## Define a PCA through the formula interface
	## Largely inspired from prcomp.formula
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
        stop("PCA applies only to numerical variables")
    na.act <- attr(mf, "na.action")
    mt <- attr(mf, "terms")
    attr(mt, "intercept") <- 0L
    x <- model.matrix(mt, mf)
    res <- pcomp.default(x, ...)
    cl[[1L]] <- as.name("pcomp")
    res$call <- cl
    if (!is.null(na.act)) {
        res$na.action <- na.act
        if (!is.null(sc <- res$x)) 
            res$x <- napredict(na.act, sc)
    }
    return(res)
}

pcomp.default <- function (x, method = c("svd", "eigen"), scores = TRUE, 
center = TRUE, scale = TRUE, tol = NULL, covmat = NULL,
subset = rep(TRUE, nrow(as.matrix(x))), ...)
{
	## Perform a PCA, either using prcomp (method = "svd"), or princomp ("eigen")
	svd.pca <- function (x, retx, center, scale, tol, ...) {
		pca <- prcomp(x, retx = retx, center = center, scale = scale, tol = tol,
			...)
		## Rework the result to make it fit in the "pcomp" object
		names(pca$sdev) <- paste("PC", 1:length(pca$sdev), sep = "") 
		if (isTRUE(!pca$center)) {
			pca$center <- rep(0, length(pca$sdev))
			names(pca$center) <- colnames(pca$rotation)
		}
		if (isTRUE(!pca$scale)) {
			pca$scale <- rep(1, length(pca$sdev))
			names(pca$scale) <- colnames(pca$rotation)
		}
		rn <- rownames(x)
			if (is.null(rn)) {
				rownames(pca$x) <- as.character(1:nrow(pca$x))
			} else {
				rownames(pca$x) <- rn
			}
		res <- list(
			loadings = structure(pca$rotation, class = "loadings"),
			scores = if (is.null(pca$x)) NULL else as.data.frame(pca$x),
			sdev = pca$sdev,
			totdev = sqrt(sum(pca$sdev^2)),
			n.obs = nrow(pca$x),
			center = pca$center,
			scale = pca$scale,
			method = "svd"
		)
		return(res)
	}
	
	eigen.pca <- function (x, cor, scores, covmat = NULL, subset, ...) {
		if (is.null(covmat)) {
			pca <- princomp(x, cor = cor, scores = scores, subset = subset, ...)
		} else {
			pca <- princomp(cor = cor, scores = scores, covmat = covmat,
				subset = subset, ...)
		}
		n <- length(pca$sdev)
		pc <- paste("PC", 1:n, sep = "")  # rename Comp.1, ... in PC1, ...
		names(pca$sdev) <- pc
		colnames(pca$loadings) <- pc
		if (!is.null(pca$scores)) {
			colnames(pca$scores) <- pc
			## If there are rownames to x, use it
			rn <- rownames(x)
			if (is.null(rn)) {
				rownames(pca$scores) <- as.character(1:nrow(pca$scores))
			} else {
				rownames(pca$scores) <- rn
			}
			pca$scores <- as.data.frame(pca$scores)
		}
		res <- list(
			loadings = pca$loadings,
			scores = pca$scores,
			sdev = pca$sdev,
			totdev = sum(pca$sdev),
			n.obs = pca$n.obs,
			center = pca$center,
			scale = pca$scale,
			method = "eigen"			
		)
		return(res)
	}

	cl <- match.call()
    cl[[1L]] <- as.name("pcomp")
	
	## Check that all variables are numeric (otherwise, issue a clear message)!
	x <- as.data.frame(x)
	if (!all(sapply(x, is.numeric)))
		stop("Cannot perform a PCA: one or more variables are not numeric.")
	
	method <- match.arg(method)
	if (method == "eigen" && !isTRUE(center))
		warning("For method 'eigen', center is always TRUE")
	res <- switch(method,
		svd = svd.pca(x, retx = scores, center = center, scale = scale,
			tol = tol, ...),
		eigen = eigen.pca(x, cor = scale, scores = scores, covmat = covmat,
			subset = subset, ...),
		stop("method must be either 'svd' or 'eigen'")
	)
	## Add a call item
	res$call <- cl
	## We return a specific object, but it is compatible (i.e., overloads), both
	## "pca" in package labdsv and "princomp" in package stats
	class(res) <- c("pcomp", "pca", "princomp")
	return(res)
}

## print method (similar to print.princomp, but reports variances instead of sds)
print.pcomp <- function (x, ...)
{
    cat("Call:\n")
    dput(x$call, control = NULL)
    cat("\nVariances:\n")
    print(x$sdev^2, ...)
	cat("\n", length(x$scale), " variables and ", x$n.obs, "observations.\n")
    invisible(x)
}

## summary method (same as summary.princomp, but with TRUE for loadings)
summary.pcomp <- function (object, loadings = TRUE, cutoff = 0.1, ...) 
{
    object$cutoff <- cutoff
    object$print.loadings <- loadings
    class(object) <- "summary.pcomp"
    object
}

## print method for summary.pcomp object (slightly modified from princomp)
print.summary.pcomp <- function (x, digits = 3, loadings = x$print.loadings,
cutoff = x$cutoff, ...) 
{
    vars <- x$sdev^2
    vars <- vars/sum(vars)
    cat("Importance of components (eigenvalues):\n")
    print(rbind(`Variance` = round(x$sdev^2, 5),
		`Proportion of Variance` = round(vars, 5), 
        `Cumulative Proportion` = round(cumsum(vars), 5)), digits = digits, ...)
    if (loadings) {
        cat("\nLoadings (eigenvectors, rotation matrix):\n")
        cx <- format(round(x$loadings, digits = digits))
        cx[abs(x$loadings) < cutoff] <- paste(rep(" ", nchar(cx[1, 1],
			type = "w")), collapse = "")
        print(cx, quote = FALSE, ...)
    }
    invisible(x)
}

## plot method
## TODO: same mechanism as for plot.lm: multiplot allowed!
plot.pcomp <- function (x,
which = c("screeplot", "loadings", "correlations", "scores"), choices = 1L:2L,
col = par("col"), bar.col = "gray", circle.col = "gray", ar.length = 0.1,
pos = NULL, labels = NULL, cex = par("cex"),
main = paste(deparse(substitute(x)), which, sep = " - "), xlab, ylab, ...)
{	
	plotScores <- function (x, choices, col, circle.col, labels, cex, main,
		xlab, ylab, ...)
	{
		if (is.null(x$scores))
			stop("no scores are available: refit with 'scores = TRUE'")
		if (is.null(labels)) {
			labels <- rownames(x$scores)
			if (is.null(labels))  # If still no labels
				labels <- as.character(1:nrow(x$scores))
		} else if (!isTRUE(!as.numeric(labels)))
			labels <- as.character(labels)
		scores <- scores(x)[, choices]
		plot(scores, type = "n", asp = 1, main = main, xlab = xlab, ylab = ylab)
		abline(h = 0, col = circle.col)
		abline(v = 0, col = circle.col)
		if (!isTRUE(!as.numeric(labels)))
			text(scores, labels = labels, col = col, cex = cex, ...)
	}
	
	which <- match.arg(which)
	main <- main[1]
	## Calculate default xlab and ylab
	labs <- paste(names(x$sdev), " (", round((x$sdev^2 / x$totdev^2) * 100,
		digits = 1), "%)", sep = "")
	if (missing(xlab)) xlab <- labs[choices[1]] else xlab
	if (missing(ylab)) ylab <- labs[choices[2]] else ylab
	switch(which,
		screeplot = screeplot(x, col = bar.col, main = main, ...),
		loadings = vectorplot(loadings(x), choices = choices, col = col,
			circle.col = circle.col, ar.length = ar.length, pos = pos, cex = cex,
			labels = if (is.null(labels)) rownames(loadings(x)) else labels,
			main = main, xlab = xlab, ylab = ylab, ...),
		correlations = vectorplot(correlation(x), choices = choices, col = col,
			circle.col = circle.col, ar.length = ar.length, pos = pos, cex = cex,
			labels = if (is.null(labels)) rownames(loadings(x)) else labels,
			main = main, xlab = xlab, ylab = ylab, ...),
		scores = plotScores(x, choices = choices, col = col, cex = cex,
			circle.col = circle.col, labels = labels, main = main,
			xlab = xlab, ylab = ylab, ...),
		stop("unknown graph type")
	)
}

## screeplot method (add cumulative variance curve to the plot)
screeplot.pcomp <- function (x, npcs = min(10, length(x$sdev)),
type = c("barplot", "lines"), col = "cornsilk", main = deparse(substitute(x)),
...) 
{
    main
    type <- match.arg(type)
    pcs <- x$sdev^2
    xp <- seq_len(npcs)
    if (type == "barplot") 
        barplot(pcs[xp], names.arg = names(pcs[xp]), main = main, 
            ylab = "Variances", col = col, ...)
    else {
        plot(xp, pcs[xp], type = "b", axes = FALSE, main = main, 
            xlab = "", ylab = "Variances", ...)
        axis(2)
        axis(1, at = xp, labels = names(pcs[xp]))
    }
    return(invisible())
}

## points method
# This is supposed to add points to a graph of scores
points.pcomp <- function (x, choices = 1L:2L, type = "p", pch = par("pch"),
col = par("col"), bg = par("bg"), cex = par("cex"), ...)
{
	if (is.null(x$scores))
		stop("no scores are available: refit with 'scores = TRUE'")
	points(scores(x)[, choices], type = type, pch = pch, col = col, bg = bg,
		cex = cex, ...)
}

## lines method
# Uses groups to draw either polygons or ellipses for each group
lines.pcomp <- function (x, choices = 1L:2L, groups, type = c("p", "e"),
col = par("col"), border = par("fg"), level = 0.9, ...)
{
	
	polygons <- function (scores, groups, n, col, border, ...) {
		for (i in 1:n) {
			sc <- na.omit(scores[as.numeric(groups) == i, ])
			if (NROW(sc) > 1) {
				pts <- chull(sc)
				## Close polygon
				pts <- c(pts, pts[1])
				polygon(sc[pts, 1], sc[pts, 2], col = col[i],
					border = border[i], ...)
			}
		}
	}

	ellipses <- function (scores, groups, n, col, border, level, ...) {
		for (i in 1:n) {
			sc <- na.omit(scores[as.numeric(groups) == i, ])
			if (NROW(sc) > 1) {
				x <- sc[, 1]
				y <- sc[, 2]
				polygon(ellipse(cor(x, y), scale = c(sd(x), sd(y)),
					centre = c(mean(x), mean(y)), level = level), col = col[i],
					border = border[i], ...)
			}
		}
	}
	
	if (is.null(x$scores))
		stop("no scores are available: refit with 'scores = TRUE'")
	if (missing(groups))
		stop("you must provide groups")
	scores <- x$scores[, choices]
	groups <- as.factor(groups)
	n <- length(levels(groups))
	col <- rep(col, length.out = n)
	border <- rep(border, length.out = n)
	type <- match.arg(type)
	switch(type,
		p = polygons(scores, groups = groups, n = n, col = col,
			border = border, ...),
		e = ellipses(scores, groups = groups, n = n, col = col,
			border = border, level = level, ...),
		stop("unknown type, currently only 'p' for polygons et 'e' for ellipses")
	)
}

## text method
text.pcomp <- function (x, choices = 1L:2L, labels = NULL, col = par("col"),
cex = par("cex"), pos = NULL, ...) {
	if (is.null(x$scores))
		stop("no scores are available: refit with 'scores = TRUE'")
	if (is.null(labels))
		labels <- as.character(1:nrow(x$scores))
	text(x$scores[, choices], labels = labels, col = col, cex = cex,
		pos = pos, ...)
}

## biplot method (note: it plots loadings, not correlations!)
biplot.pcomp <- function (x, choices = 1L:2L, scale = 1, pc.biplot = FALSE, ...) 
{
    if (length(choices) != 2) 
        stop("length of choices must be 2")
    if (!length(scores <- x$scores)) 
        stop(gettextf("object '%s' has no scores", deparse(substitute(x))), 
            domain = NA)
	if (is.complex(scores)) 
        stop("biplots are not defined for complex PCA")
    lam <- x$sdev[choices]
    n <- NROW(scores)
    lam <- lam * sqrt(n)
    if (scale < 0 || scale > 1) 
        warning("'scale' is outside [0, 1]")
    if (scale != 0) 
        lam <- lam^scale
    else lam <- 1
    if (pc.biplot) 
        lam <- lam/sqrt(n)
    biplot(t(t(scores[, choices])/lam),
		t(t(x$loadings[, choices]) * lam), ...)
    return(invisible())
}

## .panel.individuals required by pairs.pcomp
.panel.individuals <- function (x, y, ...) {
	## x and y are c(indivs, NaN, vars) => collect indivs
	pos <- 1:(which(is.nan(x))[1] - 1)
	points(x[pos], y[pos], ...)
}

## .panel.variables required by pairs.pcomp
.panel.variables <- function (x, y, ar.labels, ar.col, ar.cex, labels, col,
cex, ...)
{
	## x and y are c(indivs, NaN, vars) => collect indivs
	pos <- (which(is.nan(x))[1] + 1):length(x)
	par(new = TRUE)
	## We want to invert position of x and y here to get same one as indivs
	vectorplot(y[pos], x[pos], axes = FALSE, labels = ar.labels, col = ar.col,
		cex = ar.cex, ...)
}

## pairs plot for pcomp objects
pairs.pcomp <- function (x, choices = 1L:3L,
type = c("loadings", "correlations"), col = par("col"), circle.col = "gray",
ar.col = par("col"), ar.length = 0.05, pos = NULL, ar.cex = par("cex"),
cex = par("cex"), ...)
{
	type <- match.arg(type)
	X <- scores(x)[, choices]
	## Calculate labels
	labs <- paste(names(x$sdev), " (", round((x$sdev^2 / x$totdev^2) * 100,
		digits = 1), "%)", sep = "")[choices]
	## Add a row of NaN to separate indivs and vars
	X <- rbind(X, rep(NaN, length(choices)))
	## Add vars
	vars <- switch(type,
		loadings = loadings(x)[, choices],
		correlations = correlation(x)[, choices]
	)
	X <- rbind(X, vars)
	## Change names
	names(X) <- labs
	## Why do I get warning with non par arguments?!
	suppressWarnings(pairs(X, lower.panel = .panel.individuals,
		upper.panel = .panel.variables,
		col = col, circle.col = circle.col, ar.col = ar.col, ar.cex = ar.cex,
		ar.length = ar.length, ar.labels = rownames(vars), pos = pos,
		cex = cex, ...))	
}

## predict method
predict.pcomp <- function (object, newdata, dim = length(object$sdev), ...) 
{
    if (dim > length(object$sdev)) {
        warning("Only", length(object$sdev), " axes available\n")
        dim <- length(object$sdev)
    }
	if (missing(newdata)) 
        if (!is.null(object$scores)) 
            return(object$scores[, 1:dim])
        else stop("no scores are available: refit with 'scores = TRUE'")
    if (length(dim(newdata)) != 2L) 
        stop("'newdata' must be a matrix or data frame")
    nm <- rownames(object$loadings)
    if (!is.null(nm)) {
        if (!all(nm %in% colnames(newdata))) 
            stop("'newdata' does not have named columns matching one or more of the original columns")
        newdata <- newdata[, nm, drop = FALSE]
    } else {
        if (NCOL(newdata) != NROW(object$loadings)) 
            stop("'newdata' does not have the correct number of columns")
    }
    scale(newdata, object$center, object$scale) %*% object$loadings[, 1:dim]
}

## Extract correlation from a pcomp object. If newvar is provided, it
## calculates correlations between this new variable and corresponding PCs
## (providing that scores were calculated, and that the nrow() of new
## variable is the same as nrow(scores), assumed to be the same individuals
## as in the original PCA)
## It creates a 'corr' object
correlation.pcomp <- function (x, newvars, dim = length(x$sdev), ...)
{
	Call <- match.call()
	
	dim <- as.integer(dim)[1]
	if (dim > length(x$sdev)) {
        warning("Only", length(x$sdev), " axes available\n")
        dim <- length(x$sdev)
    }
	dims <- 1:dim
	
	if (missing(newvars)) {
		## Just extract correlations (calculated after loadings)
		if (is.null(loads <- loadings(x))) {
			return(NULL)
		} else {
			res <- sweep(loads[, dims], 2,  x$sdev[dims], "*")
			## Create a 'corr' object with this
			attr(res, "method") <- "PCA variables and components correlation"
			attr(res, "call") <- Call
			class(res) <- c("correlation", "matrix")
			return(res)
		}
	} else {
		## Calculate correlation of new variables with PCs
		## Must have same number of observations as in scores, otherwise, we got
		## the error message: "incompatible dimensions"
		if (is.null(scores <- x$scores))
			stop("no scores are available: refit with 'scores = TRUE'")
		## TODO: if these are rownames, check that they match
		res <- correlation(newvars, scores[, dims])
		## Just change method attribute
		attr(res, "method") <- "PCA variables and components correlation"
		return(res)
	}
}

## A generic function compatible with the corresponding one in labdsv package
scores <- function (x, ...)
	UseMethod("scores")

## Borrowed from scores in labdsv
## but return a data frame instead of a matrix
## TODO: check labels length, dim perhaps not the best argument name
scores.pcomp <- function (x, labels = NULL, dim = length(x$sdev), ...)
{
	if (dim > length(x$sdev)) {
        warning("Only", length(x$sdev), " axes available\n")
        dim <- length(x$sdev)
    }
	if (is.null(x$scores))
		stop("no scores are available: refit with 'scores = TRUE'")
    if (!is.null(labels)) {
		res <- as.data.frame(cbind(x$scores[, 1:dim], labels))
    } else {
        res <- as.data.frame(x$scores[, 1:dim])
    }
	return(res)
}
