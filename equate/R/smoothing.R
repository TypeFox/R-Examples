#' Frequency Distribution Presmoothing
#' 
#' These functions are used to smooth frequency distributions.
#' 
#' Loglinear smoothing is a flexible procedure for reducing irregularities in a
#' frequency distribution prior to equating, where the degree of each
#' polynomial term determines the specific moment of the observed distribution
#' that is preserved in the fitted distribution (see below for examples). The
#' \code{loglinear} function is a wrapper for \code{\link{glm}}, and is used to
#' simplify the creation of polynomial score functions and the fitting and
#' comparing of multiple loglinear models.
#' 
#' \code{scorefun}, if supplied, must contain at least one score function of
#' the scale score values. Specifying a list to \code{degrees} is an
#' alternative to supplying \code{scorefun}. Each list element in
#' \code{degrees} should be a vector equal in length to the number of variables
#' contained in \code{x}; there should also be one such vector for each
#' possible level of interaction between the variables in \code{x}.
#' 
#' For example, the default \code{degrees = list(4, 2, 2)} is recycled to
#' produce \code{list(c(4, 4, 4), c(2, 2, 2), c(2, 2, 2))}, resulting in
#' polynomials to the fourth power for each univariate distribution, to the
#' second power for each two-way interaction, and to the second power for the
#' three-way interaction.
#' 
#' Terms can also be specified with \code{grid}, which is a matrix with each
#' row containing integers specifying the powers for each variable at each
#' interaction term, including main effects. For example, the main effect to
#' the first power for the total score in a bivariate distribution would be
#' \code{c(1, 0)}; the interaction to the second power would be \code{c(2, 2)}.
#' 
#' \code{stepup} is used to run nested models based on subsets of the columns
#' in \code{scorefun}. Output will correspond to models based on columns 1 and
#' 2, 1 through 3, 1 through 4, to 1 through \code{ncol(scorefun)}. This list
#' of polynomial terms is then used to create a \code{grid} using
#' \code{expand.grid}. The \code{grid} can also be supplied directly, in which
#' case \code{degrees} will be ignored.
#' 
#' \code{compare} returns output as an \code{anova} table, comparing model fit
#' for all the models run with \code{stepup = TRUE}, or by specifying more than
#' one model in \code{models}. When \code{choose = TRUE}, the arguments
#' \code{choosemethod} and \code{chip} are used to automatically select the
#' best-fitting model based on the \code{anova} table from running
#' \code{compare}.
#' 
#' The remaining smoothing methods make adjustments to scores with low or zero
#' frequencies. \code{smoothmethod = "bump"} adds the proportion \code{jmin} to
#' each score point and then adjusts the probabilities to sum to 1.
#' \code{smoothmethod = "average"} replaces frequencies falling below the
#' minimum \code{jmin} with averages of adjacent values.
#' 
#' @param x either an object of class \dQuote{\code{freqtab}} specifying a
#' univariate or multivariate score distribution, or a \dQuote{\code{formula}}
#' object.
#' @param smoothmethod character string indicating the smoothing method to be
#' used by \code{presmoothing}. \code{"none"} returns unsmoothed frequencies,
#' \code{"bump"} adds a small frequency to each score value, \code{"average"}
#' imputes small frequencies with average values, and \code{"loglinear"} fits
#' loglinear models. See below for details.
#' @param jmin for \code{smoothmethod = "average"}, the minimum frequency, as
#' an integer, below which frequencies will be replaced (default is 1). for
#' \code{smoothmethod = "bump"}, the value to be added to each score point (as
#' a probability, with default 1e-6).
#' @param asfreqtab logical, with default \code{TRUE}, indicating whether or
#' not a frequency table should be returned. For \code{smoothmethod =
#' "average"} and \code{smoothmethod = "bump"}, the alternative is a vector of
#' frequencies. For \code{loglinear}, there are other options.
#' @param data an object of class \dQuote{\code{freqtab}}.
#' @param scorefun matrix of score functions used in loglinear presmoothing,
#' where each column includes a transformation of the score scale or
#' interactions between score scales. If missing, \code{degrees} and
#' \code{xdegree} will be used to construct polynomial score functions.
#' @param degrees list of integer vectors, each one indicating the maximum
#' polynomial score transformations to be computed for each variable at a given
#' order of interactions. Defaults (\code{degrees = list(4, 2, 2)}) are
#' provided for up to trivariate interactions. \code{degrees} are ignored if
#' \code{scorefun} or \code{grid} are provided. See below for details.
#' @param grid matrix with one column per margin in \code{x} and one row per
#' term in the model. See below for details.
#' @param rmimpossible integer vector indicating columns in \code{x} to be used
#' in removing impossible scores before smoothing, assuming internal anchor
#' variables. Impossible scores are kept by default. See below.
#' @param models integer vector indicating which model terms should be grouped
#' together when fitting multiple nested models. E.g., \code{models = c(1, 1,
#' 2, 3)} will compare three models, with the first two terms in model one, the
#' third term added in model two, and the fourth in model three.
#' @param stepup logical, with default \code{FALSE}, indicating whether or not
#' multiple nested models should be automatically fit. If \code{TRUE} and
#' \code{models} is missing, an attempt will be made to create it using
#' \code{grid} and/or \code{degrees}. Otherwise, in the absence of
#' \code{models}, each column in \code{scorefun} will define a new sequential
#' model.
#' @param compare logical, with default \code{FALSE}, indicating whether or not
#' fit for nested models should be compared. If \code{TRUE}, \code{stepup} is
#' also set to \code{TRUE} and only results from the model fit comparison are
#' returned, that is, \code{verbose} is ignored.
#' @param choose logical, with default \code{FALSE}, indicating whether or not
#' the best-fitting model should be returned after comparing fit of nested
#' models. Useful for automating model selection in simulations.
#' @param choosemethod string, with default \code{"chi"}, indicating the method
#' for selecting a best-fitting model when \code{choose = TRUE}. \code{"chi"}
#' selects the most complex model with chi-square p-value below the criterion
#' in \code{chip}. \code{"aic"} and \code{"bic"} choose the model with lowest
#' AIC and BIC, respectively. See \code{\link{anova.glm}} for details on fit
#' comparisons.
#' @param chip proportion specifying the type-I error rate for model selection
#' based on \code{choosemethod = "chi"}.
#' @param verbose logical, with default \code{FALSE}, indicating whether or not
#' full \code{glm} output should be returned.
#' @param \dots further arguments passed to other methods. For
#' \code{presmoothing}, these are passed to \code{loglinear} and include those
#' listed above.
#' @return When \code{smoothmethod = "average"} or \code{smoothmethod =
#' "bump"}, either a smoothed frequency vector or table is returned. Otherwise,
#' \code{loglinear} returns the following: \item{}{when \code{compare = TRUE},
#' an anova table for model fit} \item{}{when \code{asfreqtab = TRUE}, a
#' smoothed frequency table} \item{}{when \code{choose = TRUE}, a smoothed
#' frequency table with attribute "anova" containing the model fit table for
#' all models compared} \item{}{when \code{verbose = TRUE}, full \code{glm}
#' output, for all nested models when \code{stepup = TRUE}} \item{}{when
#' \code{stepup = TRUE} and \code{verbose = FALSE}, a \code{data.frame} of
#' fitted frequencies, with one column per model}
#' @author Anthony Albano \email{tony.d.albano@@gmail.com}
#' @seealso \code{\link{glm}}, \code{\link{loglin}}
#' @references Holland, P. W., and Thayer, D. T. (1987). \emph{Notes on the use
#' of log-linear models for fitting discrete probability distributions} (PSR
#' Technical Rep. No. 87-79; ETS RR-87-31). Princeton, NJ: ETS.
#' 
#' Holland, P. W., and Thayer, D. T. (2000). Univariate and bivariate loglinear
#' models for discrete test score distributions. \emph{Journal of Educational
#' and Behavioral Statistics, 25}, 133--183.
#' 
#' Moses, T., and Holland, P. W. (2008). \emph{Notes on a general framework for
#' observed score equating} (ETS Research Rep. No. RR-08-59). Princeton, NJ:
#' ETS.
#' 
#' Wang, T. (2009). Standard errors of equating for the percentile rank-based
#' equipercentile equating with log-linear presmoothing. \emph{Journal of
#' Educational and Behavioral Statistics, 34}, 7--23.
#' @keywords smooth models
#' @examples
#' 
#' set.seed(2010)
#' x <- round(rnorm(1000, 100, 15))
#' xscale <- 50:150
#' xtab <- freqtab(x, scales = xscale)
#' 
#' # Adjust frequencies
#' plot(xtab, y = cbind(average = freqavg(xtab),
#'   bump = freqbump(xtab)))
#' 
#' # Smooth x up to 8 degrees and choose best fitting model
#' # based on aic minimization
#' xlog1 <- loglinear(xtab, degrees = 8,
#'   choose = TRUE, choosemethod = "aic")
#' plot(xtab, as.data.frame(xlog1)[, 2],
#'   legendtext = "degree = 3")
#' 
#' # Add "teeth" and "gaps" to x
#' # Smooth with formula interface
#' teeth <- c(.5, rep(c(1, 1, 1, 1, .5), 20))
#' xttab <- as.freqtab(cbind(xscale, c(xtab) * teeth))
#' xlog2 <- presmoothing(~ poly(total, 3, raw = TRUE),
#'   xttab, showWarnings = FALSE)
#' 
#' # Smooth xt using score functions that preserve 
#' # the teeth structure (also 3 moments)
#' teeth2 <- c(1, rep(c(0, 0, 0, 0, 1), 20))
#' xt.fun <- cbind(xscale, xscale^2, xscale^3)
#' xt.fun <- cbind(xt.fun, teeth2, xt.fun * teeth2)
#' xlog3 <- loglinear(xttab, xt.fun, showWarnings = FALSE)
#' 
#' # Plot to compare teeth versus no teeth
#' op <- par(no.readonly = TRUE)
#' par(mfrow = c(3, 1))
#' plot(xttab, main = "unsmoothed", ylim = c(0, 30))
#' plot(xlog2, main = "ignoring teeth", ylim = c(0, 30))
#' plot(xlog3, main = "preserving teeth", ylim = c(0, 30))
#' par(op)
#' 
#' # Bivariate example, preserving first 3 moments of total
#' # and anchor for x and y, and the covariance
#' # between anchor and total
#' # see equated scores in Wang (2009), Table 4
#' xvtab <- freqtab(KBneat$x, scales = list(0:36, 0:12))
#' yvtab <- freqtab(KBneat$y, scales = list(0:36, 0:12))
#' Y <- as.data.frame(yvtab)[, 1]
#' V <- as.data.frame(yvtab)[, 2]
#' scorefun <- cbind(Y, Y^2, Y^3, V, V^2, V^3, V*Y)
#' wang09 <- equate(xvtab, yvtab, type = "equip",
#'   method = "chained", smooth = "loglin",
#'   scorefun = scorefun)
#' wang09$concordance
#' 
#' # Removing impossible scores has essentially no impact
#' xvlog1 <- loglinear(xvtab, scorefun, asfreqtab = FALSE)
#' xvlog2 <- loglinear(xvtab, scorefun, rmimpossible = 1:2)
#' plot(xvtab, cbind(xvlog1,
#' 	xvlog2 = as.data.frame(xvlog2)[, 3]))
#'
#' @export
presmoothing <- function(x, ...) UseMethod("presmoothing")

#' @describeIn presmoothing Default method for frequency tables.
#' @export
presmoothing.default <- function(x, smoothmethod = c("none",
	"average", "bump", "loglinear"), jmin,
	asfreqtab = TRUE, ...) {

	smoothmethod <- match.arg(smoothmethod)
	if (smoothmethod == "none")
		return(x)
	else if (smoothmethod == "average")
		return(freqavg(x, jmin = jmin,
			asfreqtab = asfreqtab))
	else if (smoothmethod == "bump")
		return(freqbump(x, jmin = jmin,
			asfreqtab = asfreqtab))
	else if (smoothmethod == "loglinear")
		return(loglinear(x, asfreqtab = asfreqtab, ...))
}

#----------------------------------------------------------------
# Formula method

#' @describeIn presmoothing Method for \dQuote{\code{formula}} objects.
#' @export
presmoothing.formula <- function(x, data, ...) {
	
	formula <- terms(as.formula(x))
	if (attr(formula, "response") == 0) {
		formula <- terms(reformulate(attr(formula,
			"term.labels"), "count"))
	}
	if (attr(formula, "intercept") == 0) {
		attributes(formula)$intercept <- 1
		warning("an intercept has been added to the model")
	}

	scorefun <- as.data.frame(model.matrix(formula,
		data = as.data.frame(data))[, -1])
	loglinear(data, scorefun, ...)
}


#----------------------------------------------------------------
# Exported internal function for loglinear smoothing

#' @rdname presmoothing
#' @export
loglinear <- function(x, scorefun, degrees = list(4, 2, 2), grid,
	rmimpossible, asfreqtab = TRUE, models,
	stepup = !missing(models), compare = FALSE, choose = FALSE,
	choosemethod = c("chi", "aic", "bic"), chip,
	verbose = FALSE, ...) {

	# Powers in higher order interactions should never 
	# be larger than in lower - they will be ignored
	
	xd <- as.data.frame(x)
	nx <- ncol(xd) - 1
	
	# Remove impossible scores, assuming internal anchors
	if (!missing(rmimpossible) && is.numeric(rmimpossible)) {
		keepi <- apply(xd[, sort(rmimpossible)], 1,
			function(y) all(y[-1] <= y[1]))
		xd <- xd[keepi, ]
	}
	else
		keepi <- rep(TRUE, nrow(xd))
	if (choose) compare <- TRUE
	if (missing(scorefun))
		scorefun <- sf(xd, degrees, grid, stepup, compare)
		if (stepup | compare) {
			models <- attributes(scorefun)$models
			mnames <- attributes(scorefun)$mnames
		}
	else {
		scorefun <- scorefun[keepi, ]
		if (stepup | compare) {
			if (missing(models))
				models <- 1:ncol(scorefun)
			mnames <- unique(models)
		}
	}
	if (nrow(scorefun) != nrow(xd))
		stop("'scorefun' must contain the same ",
			"number of rows as 'x'")
	scorefun <- data.frame(f = xd[, nx + 1],
		scorefun, check.names = FALSE)
	if (stepup | compare) {
		if (ncol(scorefun) < 3)
			stop(paste("cannot run multiple models with only",
				ncol(scorefun) - 1, "model terms"))
		snames <- colnames(scorefun)[-1]
		out <- lapply(unique(models), function(i)
			glm(scorefun[, c("f", snames[models <= i])],
				family = poisson))
		names(out) <- mnames
	}
	else
		out <- glm(scorefun, family = poisson)
	if (compare) {
		nm <- length(out)
		resdf <- as.numeric(lapply(out, function(y) y$df.residual))
		resdev <- as.numeric(lapply(out, function(y) y$deviance))
		aic <- as.numeric(lapply(out, AIC))
		bic <- as.numeric(lapply(out, BIC))
		tab <- data.frame(resdf, resdev, aic, bic,
			c(NA, -diff(resdf)), c(NA, 
			-diff(resdev)))
		vars <- lapply(out, function(y) paste(deparse(formula(y)), 
			collapse = "\n"))
		dimnames(tab) <- list(1:nm, c("Resid. Df", "Resid. Dev", 
			"AIC", "BIC", "Df", "Deviance"))
		tab <- stat.anova(tab, test = "Chisq", scale = 1, 
			df.scale = Inf,
			n = length(out[[order(resdf)[1]]]$residuals))
		atab <- structure(tab,
			heading = c("Analysis of Deviance Table\n",
			paste("Model ", format(1:nm), ": ", vars,
			sep = "", collapse = "\n")),
			class = c("anova", "data.frame"))
		if (choose) {
			glmi <- glmselect(atab, choosemethod, chip)
			stab <- as.freqtab(cbind(xd[, 1:nx],
				out[[glmi]]$fitted),
				scales = scales(x, 1:nx))
			attr(stab, "anova") <- atab
			return(stab)
		} else return(atab)
	} else if (verbose)
		return(out)
	else if (stepup)
		return(data.frame(lapply(out, fitted),
			check.names = FALSE))
	else if (asfreqtab)
		return(as.freqtab(cbind(xd[, 1:nx],
			out$fitted), scales = scales(x, 1:nx)))
	else return(out$fitted)
}

#----------------------------------------------------------------
# Internal function for selecting model from anova data.frame

glmselect <- function(x, choosemethod = c("chi", "aic", "bic"),
	chip) {
	
	if (class(x)[1] != "anova")
		stop("'x' must be an anova table, output from 'glm'")
	nm <- nrow(x)
	choosemethod <- match.arg(choosemethod)
	if (choosemethod == "chi") {
		if (missing(chip))
			chip <-  1 - (1 - .05)^(1/(nm - 1))
		chib <- x[, 7] < chip
		out <- ifelse(any(chib, na.rm = T),
			max(which(chib)), 1)
	}
	else if (choosemethod == "aic")
		out <- which.min(x$AIC)
	else if (choosemethod == "bic")
		out <- which.min(x$BIC)
	return(out)
}

#----------------------------------------------------------------
# Internal function for creating score function and models

sf <- function(x, degrees, grid, stepup = FALSE, compare = stepup) {
	x <- as.data.frame(x)
	nx <- ncol(x) - 1
	if (missing(grid)) {
		if (length(degrees) < nx) # must be at least 0 for higher orders
			degrees[(length(degrees) + 1):nx] <- 0
		degrees <- lapply(degrees, function(y)
			rep(y, nx)[1:nx])
		# Start grid without intercept
		grid <- cbind(expand.grid(lapply(degrees[[1]],
			function(y) 0:y))[-1, ])
		# Remove higher order interactions as necessary
		if (nx > 1) {
			for(i in 2:nx) {
				# Make sure higher orders don't contain larger powers
				# They're already excluded in grid
				#degrees[[i]] <- pmin(degrees[[i - 1]],
				#	degrees[[i]])
				rm1 <- apply(grid, 1, function(y)
					sum(y == 0) == (nx - i))
				rm2 <- apply(sapply(1:nx, function(j)
					grid[, j] > degrees[[i]][j]), 1, any)
				grid <- grid[!(rm1 & rm2), ]
			}
			# Sort grid by orders
			grid <- cbind(grid[order(apply(grid, 1, function(y)
				sum(y == 0)), decreasing = T), ])
			os <- factor(nx - apply(grid, 1, function(y)
				sum(y == 0)))
			grid <- do.call("rbind", by(grid, os,
				function(y) y[order(apply(y, 1, max)), ]))
		}
		else
			os <- factor(nx - apply(grid, 1, function(y)
				sum(y == 0)))
	}
	scorefun <- NULL
	for(j in 1:nrow(grid)) {
		tempfun <- sapply(1:nx, function(k)
			x[, k]^grid[j, k])
		scorefun <- cbind(scorefun,
			apply(tempfun, 1, prod))
	}
	dimnames(scorefun) <- NULL
	colnames(scorefun) <- apply(grid, 1, paste,
		collapse = ".")
	if (stepup | compare) {
		# Create model index
		attr(scorefun, "models") <- as.numeric(factor(paste(os,
			apply(grid, 1, max), sep = ".")))
		attr(scorefun, "mnames") <- unique(paste(os,
			apply(grid, 1, max), sep = "."))
	}
	return(scorefun)
}

#----------------------------------------------------------------
# Frequency adjustment
# Bump frequencies upward by a small amount

#' @rdname presmoothing
#' @export
freqbump <- function(x, jmin = 1e-6, asfreqtab = FALSE, ...) {

	x <- as.data.frame(x)
	nc <- ncol(x)
	f <- x[, nc]/sum(x[, nc])
	out <- (f + jmin)/(1 + (max(x[, 1]) + 1)*jmin)
	out <- ((f + jmin)/sum(f + jmin))*sum(x[, nc])
	
	if (asfreqtab)
		return(as.freqtab(cbind(x[, -nc], out)))
	else
		return(out)
}

#----------------------------------------------------------------
# Frequency averaging

#' @rdname presmoothing
#' @export
freqavg <- function(x, jmin = 1, asfreqtab = FALSE, ...) {
	
	xtab <- x <- as.data.frame(x)
	nc <- ncol(xtab)
	if (nc > 2)
		stop("frequency averaging only supported ",
			"for univariate 'x'")
		
	x <- cbind(x, 0, 0, 0, 0, 0)
	ks <- 1
	while(ks <= nrow(x) & x[ks, 2] < jmin)
		ks <- ks + 1
	x[1:ks, 3] <- 1
	x[1:ks, 4] <- ks

	lls <- nrow(x)
	while(lls >= 0 & x[lls, 2] < jmin)
		lls <- lls - 1
	x[lls:nrow(x), 3] <- lls
	x[lls:nrow(x), 4] <- nrow(x)

	ss <- ks + 1
	tts <- lls - 1
	for(j in ss:tts) {
		if (x[j, 2] < jmin) {
			lls <- j
			ks <- j
			while(lls >= 1 & x[lls, 2] < jmin)
				lls <- lls - 1
			while(ks <= nrow(x) & x[ks, 2] < jmin)
				ks <- ks + 1
			x[lls:ks, 3] <- lls
			x[lls:ks, 4] <- ks
		}
	}

	for(p in 1:(nrow(x) - 1)) {
		if (x[p, 4] > 0 & x[p, 4] == x[p + 1, 3]) {
			if (x[p, 3] > 0)
				lls <- x[p, 3]
			if (x[p + 1, 4] > 0)
				ks <- x[p + 1, 4]
			x[lls:ks, 3] <- lls
			x[lls:ks, 4] <- ks
		}
	}

	for(j in 1:nrow(x)) {
		lls <- x[j, 3]
		if (lls == 0)
			lls <- j
		ks <- x[j, 4]
		if (ks == 0)
			ks <- j
		sumit <- 0
		sumit <- sumit + sum(x[lls:ks, 2])
		for(i in lls:ks) {
			x[i, 5] <- sumit
			x[i, 6] <- x[j, 4] - x[j, 3] + 1
			x[i, 7] <- x[i, 5]/x[i, 6]
		}
		j <- j + x[j, 4] - x[j, 3]
	}

	colnames(x)[c(1, 2, 6, 7)] <-
		c("score", "count", "b", "acount")

	if (asfreqtab)
		return(as.freqtab(cbind(xtab[, -nc],
			x[, 7])))
	else
		return(x[, 7])
}
