#' Observed Score Linking and Equating
#' 
#' This function links the scale of \code{x} to the scale of \code{y} for the
#' single-group, equivalent-groups, and nonequivalent-groups with anchor test
#' designs. A summary method is also provided.
#' 
#' Equating is typically performed on two frequency tables, \code{x} and
#' \code{y}. In this case, the scores from both are used to define the equating
#' function that links the scale of \code{x} to \code{y}. The equivalent-groups
#' design is assumed when \code{x} and \code{y} are separate, univariate
#' frequency tables. The nonequivalent-groups design is assumed when a
#' \code{method} is specified, and \code{x} and \code{y} are separate
#' multivariate frequency tables. Finally, a single-group design is assumed
#' when \code{x} is a bivariate frequency table (containing scores on X and Y)
#' and \code{y} is missing.
#' 
#' The single-group design currently only differs from the equivalent groups
#' design in that presmoothing can be used to preserve bivariate moments for
#' \code{x} and \code{y} in the single-group design, whereas in the
#' equivalent-groups design, with \code{x} and \code{y} specified separately,
#' presmoothing is performed separately. If presmoothing is not performed via
#' \code{equate}, the single-group and equivalent-groups designs produce the
#' same result.
#' 
#' When \code{x} is a vector of scores and equating output is supplied for
#' \code{y}, no other arguments are required. Scores from \code{x} are
#' converted directly to the scale indicated in \code{y}. If \code{y} is a
#' composite equating, composite equated scores will be returned based on the
#' weighted combination of equating functions included in \code{y}.
#' 
#' When \code{x} is a \code{list} of frequency tables, each element in
#' \code{args} must be a named list of equating arguments. In this case, the
#' length of \code{args} corresponds to the number of equatings that will be
#' performed. The arguments for each equating are specified as they would be
#' when \code{x} and \code{y} are frequency tables, except for \code{x} and
#' \code{y}; the frequency tables to be equated are specified in \code{args} by
#' referencing their names in the list of frequency tables. See below for
#' examples.
#' 
#' Six equating types are currently supported: identity, mean, linear, and
#' equipercentile, as described by Kolen and Brennan (2004); circle-arc
#' equating, as described by Livingston and Kim (2009); and a general linear
#' function that extends the traditional identity, mean, and linear types.
#' Corresponding linking methods are also supported. The equating design is
#' implied by the \code{method} argument, where \code{"none"} (default)
#' indicates that no method is needed (because examinees taking forms X and Y
#' are assumed to be the same or equivalent). The nominal weights, Tucker,
#' Levine observed score, Levine true score, frequency estimation,
#' Braun/Holland, and chained equating methods are supported for the
#' nonequivalent-groups with anchor test design. All but the Levine true score
#' and chained method rely on a \dQuote{synthetic} distribution of scores
#' (Braun and Holland, 1982), a weighted combination of \code{x} and \code{y}.
#' 
#' Depending on the equating method, the following additional arguments may be
#' required:
#' \describe{ \item{midp}{ coordinates for the midpoint of
#' the equating line, used in general linear and circle-arc equating. }
#' \item{cx, cy, sx, sy}{ parameters used in general linear equating.
#' See below for details. } \item{wax, way, wbx, wby}{ weights used
#' when finding the slope and intercept in general linear equating. See below.
#' } \item{ws}{ value between 0 and 1 specifying the weight applied to
#' form X scores (and implicitly specifying the form Y weight as \code{1 - ws})
#' when estimating the synthetic population. When set to -1 (the default),
#' proportional weights are calculated for X and Y based on sample size. }
#' \item{internal}{ logical indicating whether or not the anchor item
#' scores are included in the total scores. This applies only to the Levine
#' method, as all other methods assume an internal anchor test. Default is
#' \code{TRUE}. } \item{lts}{ logical indicating whether or not to use
#' levine true score (\dQuote{lts}) equating. Default is \code{FALSE}. }
#' \item{smoothmethod}{ string indicating one of four smoothing methods
#' to be used in equipercentile equating: \code{"none"} (default),
#' \code{"average"}, \code{"bump"}, and \code{"loglinear"} (see below). }
#' \item{chainmidp}{ string specifying the type of chained linear
#' equating used to obtain the midpoint in chained circle-arc equating, whether
#' \code{"mean"} (default) or \code{"linear"}. } \item{simple}{
#' logical, with default \code{TRUE}, indicating whether or not simplified
#' circle-arc equating should be used (see below). } \item{reps}{ the
#' number of replications to use in bootstrapping. Passed to
#' \code{\link{bootstrap}}. } \item{xp, yp}{ optional parametric
#' distributions, as frequency tables, replacing \code{x} and \code{y} when
#' bootstrapping. } \item{xn, yn}{ sample sizes to be drawn from
#' \code{x} and \code{y}, or \code{xp} and \code{yp}, at each bootstrap
#' replication. These default to the observed sample sizes. }
#' \item{crit}{ a vector of equated scores serving as the criterion
#' equating function when calculating bootstrap bias and RMSE; both are
#' returned when \code{crit} is specified. } }
#' 
#' General linear equating is a new
#' approach to estimating a linear linking or equating function. The slope and
#' intercept of the line are estimated based on multiple sources of
#' information, including the means and standard deviations of X and Y, and
#' other values supplied through \code{cx} and \code{cy}, representing the
#' centrality of X and Y, and \code{sx} and \code{sy}, representing the scaling
#' or variability of X and Y. The weights \code{wax} and \code{way} specify the
#' proportional weighting given to the standard deviations of X and Y, and
#' indirectly the weighting given to \code{sx} and \code{sy}, in determining
#' the slope. \code{wbx} and \code{wby} specify the proportional weighting
#' given to the means of X and Y, and indirectly the weighting given to
#' \code{cx} and \code{cy}, in determining the intercept. Synthetic means and
#' standard deviations will be used when appropriate. Chained general linear
#' equating is not currently supported.
#' 
#' For equipercentile equating under the random groups design, three smoothing
#' options are available: \code{smoothmethod = "average"} and
#' \code{smoothmethod = "bump"} require the additional argument \code{jmin},
#' and loglinear smoothing (\code{smoothmethod = "loglinear"}) requires either
#' a score function or maximum polynomial terms. For frequency estimation and
#' chained methods, only smoothing methods \code{"bump"} and \code{"loglinear"}
#' are supported. See the \code{\link{presmoothing}} function for details and
#' examples.
#' 
#' In equipercentile equating, the high point for \code{y}, i.e.,
#' \code{highp[2]}, is used to obtain form Y equivalents of form X scores with
#' percentile ranks of 100. Typically this is set to be the number of score
#' points in the form Y scale, which assumes that scores are integers ranging
#' from 1 (or 0) to the total number of items, and that each item is scored
#' correct/incorrect. Scores on other scales (such as scales which include
#' negative values, or which exclude zero) may also be used. In such cases
#' \code{highp[2]} can be set to the highest possible score on form Y, or
#' alternatively the highest observed score on Y.
#' 
#' \code{lowp} and \code{highp} are used to define the slope and intercept of
#' the identity linking function. When the score scales for X and Y are
#' equivalent, the identity function is simply the unequated X scale; however,
#' when forms differ in their scales, e.g., because of changes in content or
#' length, the identity linking function will map X onto Y based on the low and
#' high coordinates.
#' 
#' The simplified approach to circle-arc equating, as demonstrated by
#' Livingston and Kim (2009), involves combining a circle-arc with the identity
#' function. When the low and high scores differ for the X and Y scales, this
#' becomes the identity linking function. The linear component can be omitted,
#' and symmetric circle-arc equating used, with \code{simple = FALSE}. The
#' result is an equating function based only on the circle-arc that passes
#' through the points \code{lowp}, \code{highp}, and the estimated midpoint.
#' 
#' Analytical standard errors are not currently returned. With \code{boot =
#' TRUE}, bootstrap standard errors are estimated using a default of \code{reps
#' = 100} replications, sampling the maximum amount from each score
#' distribution (controlled by the arguments \code{xn} and \code{yn}). See
#' \code{\link{bootstrap}} for details and examples, including how to obtain
#' bootstrap bias and RMSE.
#' 
#' @param x,y for the default method, \code{x} must be a vector of scores and
#' \code{y} an object of class \dQuote{\code{equate}}, the output of a previous
#' equating. The standard usage is to provide \code{x} as a frequency table of
#' class \dQuote{\code{freqtab}}, where \code{y} is also a frequeny table, and
#' \code{x} is equated to \code{y}; if \code{y} is missing, a single-group
#' design is assumed. Finally, \code{x} may be a list of two or more frequency
#' tables, in which case the required arguments for one or more equatings are
#' listed in \code{args}. See below for details.
#' @param type the type of equating. See below for details.
#' @param method the equating method, where \code{"none"} (default) indicates
#' equating under the single-group or equivalent-groups design, and
#' \code{"nominal weights"}, \code{"tucker"}, \code{"levine"}, \code{"frequency
#' estimation"}, \code{"braun/holland"}, and \code{"chained"} indicate the
#' corresponding methods under the nonequivalent groups design.
#' @param name an optional name, used to label the output. If missing, a name
#' will be created using \code{type} and \code{method}.
#' @param lowp,highp two vectors, each of length 2, specifying the coordinates
#' for the low and high points of the X and Y scales. \code{lowp} defaults to
#' the minimums and \code{highp} the maximums of the scales. Recycled if
#' necessary. When \code{lowp = "obs"} or \code{highp = "obs"}, minimum and
#' maximum observed scores are used.
#' @param boot logical indicating whether or not bootstrapping should be
#' performed. Default is \code{FALSE}. See below and the
#' \code{\link{bootstrap}} function for details.
#' @param verbose logical, with default \code{TRUE}, indicating whether or not
#' full output should be returned. When \code{FALSE}, only the equated scores
#' are returned.
#' @param args list of arguments passed to \code{equate}. See below for
#' details.
#' @param object output from either an equating or list of equatings, produced
#' by the \code{equate} function.
#' @param \dots further arguments passed to or from other functions, including
#' arguments specific to different equating methods. See below for details.
#' @return When \code{y} contains output from an equating, a vector of equated
#' scores is returned. Otherwise, an object of class \dQuote{\code{equate}} is
#' returned, listing the following components, some of which are dependent on
#' the equating type, method, and smoothing: \item{name}{\code{name} for the
#' equating} \item{type}{equating type} \item{method}{equating method}
#' \item{design}{equating design, inferred from the method} \item{x,
#' y}{original frequency tables for X and Y} \item{concordance}{conversion
#' table containing scores on X with their form Y equivalents.}
#' \item{points}{low and high points defining the identity line, and midpoints
#' for general linear and circle-arc equating} \item{weight}{weights used in
#' general linear equating} \item{internal, lts, jmin, degree, xdegree,
#' scorefun}{additional arguments, as supplied in \dots{}}
#' \item{coefficients}{conversion coefficients intercept and slope; for
#' circle-arc equating, circle center points and radius are also included; for
#' general linear equating, slope and intercept components are included}
#' \item{ws}{weight applied to X in synthetic estimation}
#' \item{synthstats}{means and standard deviations for the synthetic
#' distributions} \item{xsynthetic, ysynthetic}{frequency tables for the
#' synthetic distributions} \item{smoothmethod}{smoothing method}
#' \item{xsmooth, ysmooth}{smoothed frequency tables for X and Y}
#' \item{bootstraps}{list of bootstrap standard errors, and, optionally, bias
#' and RMSE}
#' The summary method returns
#' @author Anthony Albano \email{tony.d.albano@@gmail.com}
#' @seealso \code{\link{freqbump}}, \code{\link{freqavg}},
#' \code{\link{loglinear}}, \code{\link{bootstrap}}
#' @references Kolen, M. J., and Brennan, R. L. (2004). \emph{Test Equating,
#' Scaling, and Linking}. (2nd ed.), New York: Springer.
#' 
#' Livingston, S. A., and Kim, S. (2009). The circle-arc method for equating in
#' small samples, \emph{Journal of Educational Measurement, 46}, 330--343.
#' @keywords methods
#' @examples
#' 
#' # See vignette("equatevignette") for a description of methods
#' # and additional examples
#' 
#' # Random groups equating for (1) identity, (2) mean, 
#' # (3) linear, (4) equipercentile with loglinear
#' # smoothing, and (5) a composite of mean and identity
#' rx <- as.freqtab(ACTmath[, 1:2])
#' ry <- as.freqtab(ACTmath[, c(1, 3)])
#' set.seed(2007)
#' 
#' req1 <- equate(rx, ry, type = "i", boot = TRUE, reps = 5)
#' req2 <- equate(rx, ry, type = "m", boot = TRUE, reps = 5)
#' req3 <- equate(rx, ry, type = "l", boot = TRUE, reps = 5)
#' req4 <- equate(rx, ry, type = "e", boot = TRUE, reps = 5,
#'   smooth = "loglin", degree = 3)
#' req5 <- composite(list(req1, req2), wc = .5, symmetric = TRUE)
#' 
#' # Compare equating functions
#' plot(req1, req2, req3, req4, req5[[1]], addident = FALSE)
#' 
#' # Compare boostrap standard errors
#' # Bootstrapping isn't supported for composite equating
#' plot(req1, req2, req3, req4, addident = FALSE, out = "se",
#'   legendplace = "topleft")
#' 
#' # Nonequivalent groups design for (1) Tucker linear,
#' # (2) frequency estimation , and (3) Braun/Holland linear
#' nx <- freqtab(KBneat$x, scales = list(0:36, 0:12))
#' ny <- freqtab(KBneat$y, scales = list(0:36, 0:12))
#' 
#' neq1 <- equate(nx, ny, type = "linear", method = "tuck", ws = 1)
#' neq2 <- equate(nx, ny, type = "equip", method = "freq", ws = 1)
#' neq3 <- equate(nx, ny, type = "linear", method = "braun", ws = 1)
#' 
#' # Compare equated scores
#' round(cbind(xscale = 0:36, tucker = neq1$conc$yx,
#' 	fe = neq2$conc$yx, braun = neq3$conc$yx), 2)
#' 
#' # Multiple linkings using PISA reading booklet 6
#' # clusters 3a, 5, 6, and 7
#' r3 <- freqtab(PISA$totals$b6$r3a, scales = 0:15)
#' r5 <- freqtab(PISA$totals$b6$r5, scales = 0:15)
#' r6 <- freqtab(PISA$totals$b6$r6, scales = 0:15)
#' r7 <- freqtab(PISA$totals$b6$r7, scales = 0:14)
#' eqargs <- list(r3r5 = list(type = "linear", x = "r3", y = "r5",
#'     name = "Linear Linking PISA r3 to r5"),
#'   r5r6 = list(type = "linear", x = "r5", y = "r6",
#'     name = "Linear Linking PISA r5 to r6"),
#'   r6r7 = list(type = "linear", x = "r6", y = "r7",
#'     name = "Linear Linking PISA r5 to r7"))
#' req <- equate(list(r3 = r3, r5 = r5, r6 = r6, r7 = r7), eqargs)
#' 
#' # Put PISA r3 on the scale of r7 using the linking chain
#' # Compare to a direct linking of r3 to r7
#' equate(equate(req$r3r5$conc$yx, req$r5r6), req$r6r7)
#' equate(r3, r7, "linear")$conc$yx
#' 
#' # Linking PISA cluster r3a to r5 with multiple anchors
#' m367 <- freqtab(PISA$totals$b6[1:198, c("r3a", "r6", "r7")],
#' 	scales = list(0:15, 0:16, 0:14))
#' m567 <- freqtab(PISA$totals$b6[199:396, c("r5", "r6", "r7")],
#' 	scales = list(0:15, 0:16, 0:14))
#' meq1 <- equate(m367, m567, type = "mean", method = "nom")
#' meq2 <- equate(m367, m567, type = "mean", method = "tuck")
#' meq3 <- equate(m367, m567, type = "lin", method = "tuck")
#' meq4 <- equate(m367, m567, type = "equip", method = "freq",
#' 	smooth = "log", show = FALSE)
#' meq <- equate(m367, m567, type = "mean", method = "nom")
#' plot(meq1, meq2, meq3, meq4, meq, req[[1]])
#' 
#' @export equate
equate <- function(x, ...) UseMethod("equate")

#' @describeIn equate Equating a list of frequency tables.
#' @export
equate.list <- function(x, args, ...) {
	
	if(length(x) < 2)
		stop("'x' must contain at least 2 data sets")
	if(length(x) == 2 & (is.null(names(x)) |
		any(names(x) == "")))
		names(x) <- c("x", "y")
	if(is.null(names(x)) | any(names(x) == ""))
		stop("all elements of 'x' must be named")
	if(!is.list(args))
		stop("'args' must be a list")
	if(!all(sapply(args, is.list)))
		stop("each element in 'args' must be a list")
	if(length(x) > 2) {
		if(!all(unlist(lapply(args, function(z)
			"x" %in% names(z)))))
			stop("'args' must include 'x' for each equating")
	}
	dots <- list(...)
	neq <- length(args)
	
	out <- vector("list", length = neq)
	for(i in 1:neq) {
		eqargs <- args[[i]]
		eqargs[names(dots)] <- dots
		if(is.null(eqargs$x))
			eqargs$x <- names(x)[1]
		if(is.null(eqargs$y))
			eqargs$y <- names(x)[2]
		eqargs$x <- x[[eqargs$x]]
		eqargs$y <- x[[eqargs$y]]
		out[[i]] <- do.call(equate.freqtab, eqargs)
	}
	names(out) <- names(args)

	return(as.equate.list(out))
}

as.equate.list <- function(x) {
	
	class(x) <- "equate.list"
	return(x)
}

is.equate.list <- function(x) {
	
	return(class(x)[1] == "equate.list")
}

assign("[.equate.list", function(x, i){

	out <- as.equate.list(NextMethod("["))
	return(out)
})

#' @describeIn equate Equating frequency distributions in \code{x} and \code{y}.
#' @export
equate.freqtab <- function(x, y, type = c("identity",
	"mean", "linear", "general linear", "circle-arc",
	"equipercentile"), method = c("none", "nominal weights",
	"tucker", "levine", "frequency estimation", "chained",
	"braun/holland"), name, lowp, highp, boot = FALSE,
	verbose = TRUE, ...) {
	
	if(missing(y)) 
        y <- margin(x, 2)
	if(missing(lowp))
		lowp <- c(min(scales(x)), min(scales(y)))
	else if(pmatch(lowp[1], "obs", 0))
		lowp <- c(min(x), min(y))
	else if(length(lowp) == 1)
		lowp[2] <- lowp[1]
	if(missing(highp))
		highp <- c(max(scales(x)), max(scales(y)))
	else if(pmatch(highp[1], "obs", 0))
		highp <- c(max(x), max(y))
	else if(length(highp) == 1)
		highp[2] <- highp[1]
	
	type <- match.arg(type)
	method <- match.arg(method)
	if(type %in% c("identity", "mean", "linear", "general linear"))
		eqout <- linear(x, y, type = type, method = method,
			lowp = lowp, highp = highp, verbose = verbose, ...)
	else if(type == "equipercentile")
		eqout <- equipercentile(x, y, type = type, method = method,
			lowp = lowp, highp = highp, verbose = verbose, ...)
	else if(type == "circle-arc")
		eqout <- circlearc(x, y, type = type, method = method,
			lowp = lowp, highp = highp, verbose = verbose, ...)

	if(verbose) {
		if(missing(name)) {
			name <- ifelse(method == "none", type,
				paste(method, type))
			xname <- ifelse(exists(deparse(substitute(x)), 1,
				inherits = FALSE), deparse(substitute(x)), "x")
			if(margins(y) == margins(x))
				yname <- ifelse(exists(deparse(substitute(y)), 1,
					inherits = FALSE), deparse(substitute(y)), "y")
			else yname <- xname
			name <- paste(gsub("\\b(\\w)", "\\U\\1", name,
				perl = TRUE), "Equating:", xname, "to", yname)
		}
		des <- switch(design(x), "sg" = "single group",
			"cb" = "counterbalanced",
			"eg" = "equivalent groups",
			"ng" = "nonequivalent groups")
		out <- c(list(name = name, type = type,
			method = method, design = des), eqout)
		out <- as.equate(out)
		if(boot)
			out$bootstraps <- bootstrap.equate(out, ...)
	}
	else out <- eqout

	return(out)
}

#' @describeIn equate Default equating method for a vector of raw scores \code{x}
#' and equating output in \code{y}.
#' @export
equate.default <- function(x, y, ...) {

	if(is.equate(y))
		yx <- convert(x, y, ...)
	else if(is.equate.list(y)) {
		if(is.composite(y[[1]])) {
			yx <- sapply(y[-1], function(z) convert(x, z))
			wcs <- if(!is.null(y[[1]]$wcs)) y[[1]]$wcs
				else y[[1]]$wc
			yx <- c(yx %*% wcs)
		}
		else yx <- sapply(y, function(z) convert(x, z))
	}
	else
		stop("'y' must be an 'equate' or 'equate.list' object")
		
	return(yx)
}

convert <- function(x, y, ...) {
	
	if(!is.equate(y))
		stop("'y' must be an 'equate' object")
		
	if(y$type == "equipercentile") {
		if(y$method == "none") {
			if(y$smoothmethod == "none") {
				xtab <- y$x
				ytab <- y$y
			}
			else {
				xtab <- y$xsmooth
				ytab <- y$ysmooth
			}
			p <- px(x, xtab)
		}
		else if(y$method == "frequency estimation") {
			xtab <- margin(y$xsynthetic)
			ytab <- margin(y$ysynthetic)
			p <- px(x, xtab)
		}
		else {
			if(y$smoothmethod == "none") {
				xtab <- margin(y$x)
				xvtab <- margin(y$x, 2)
				ytab <- margin(y$y)
				yvtab <- margin(y$y, 2)
			}
			else {
				xtab <- margin(y$xsmooth)
				xvtab <- margin(y$xsmooth, 2)
				ytab <- margin(y$ysmooth)
				yvtab <- margin(y$ysmooth, 2)
			}
			vx <- equip(px(x, xtab), xvtab)$yx
			p <- px(vx, yvtab)
		}
		yx <- equip(p, ytab)$yx
	}
	else if(y$type == "circle-arc") {
		yx <- circ(x, y$points[1, ], y$points[2, ],
			y$points[3, ], y$coef[1], y$coef[2], y$coef[3:4],
			y$coef[5], simple = y$simple)
	}
	else {
		yx <- lin(x, y$coef[1], y$coef[2])
	}
	names(yx) <- NULL

	return(yx)
}

as.equate <- function(x) {
	
	class(x) <- "equate"
	return(x)
}

is.equate <- function(x) {
	
	return(class(x)[1] == "equate")
}

#' @export
print.equate <- function(x, ...) {

	cat("\n")
	cat(x$name, "\n\n")
	cat("Design:", x$design, "\n\n")
	if(x$type == "equipercentile") {
		cat("Smoothing Method: ")
		cat(switch(x$smoothmethod,
			bump = "adjusted frequency presmoothing",
			average = "average frequency presmoothing",
			loglinear = "loglinear presmoothing", "none"),
			"\n\n")
	}
	if(!is.null(x$ws))
		cat("Synthetic Weighting for x:",
			x$ws, "\n\n")
	xm <- margin(x$x)
	stats <- rbind(x = summary(xm),
		y = summary(margin(x$y)),
		yx = summary(as.freqtab(cbind(x$concordance[, 2],
			c(xm)))))
	cat("Summary Statistics:\n")
		print(round(stats, 2))
		cat("\n")

	if(!is.null(x$coef)) {
		cat("Coefficients:\n")
		print(round(x$coef, 4))
		cat("\n")
	}

	invisible(x)
}

# @describeIn equate Summary method for \dQuote{"equate"} objects.
#' @rdname equate
#' @export
summary.equate <- function(object, ...) {
  
  xtab <- as.data.frame(margin(object$x))
  ytab <- as.data.frame(margin(object$y))
  if(object$type == "equipercentile" &&
      object$smoothmethod != "none") {
    xtab <- data.frame(xtab,
      smooth = c(margin(object$xsmooth)))
    ytab <- data.frame(ytab,
      smooth = c(margin(object$ysmooth)))
  }
  xvtab <- yvtab <- NULL
  if(object$design == "nonequivalent groups") {
    xvtab <- as.data.frame(margin(object$x, 2))
    yvtab <- as.data.frame(margin(object$y, 2))
    if(object$type == "equipercentile" &&
        object$smoothmethod != "none") {
      xvtab <- data.frame(xvtab,
        smooth = c(margin(object$xsmooth, 2)))
      yvtab <- data.frame(yvtab,
        smooth = c(margin(object$ysmooth, 2)))
    }
    if(object$type != "composite" && (object$method ==
        "frequency estimation" | object$method ==
        "braun/holland")) {
      xtab <- data.frame(xtab,
        synthetic = c(margin(object$xsynthetic)))
      ytab <- data.frame(ytab,
        synthetic = c(margin(object$ysynthetic)))
      xvtab <- data.frame(xvtab,
        synthetic = c(margin(object$xsynthetic, 2)))
      yvtab <- data.frame(yvtab,
        synthetic = c(margin(object$ysynthetic, 2)))
    }
  }
  yxtab <- as.freqtab(cbind(object$concordance[, 2], xtab$count))
  yxstats <- summary(yxtab)
  rownames(yxstats) <- "observed"
  
  xstats <- ystats <- NULL
  for(i in 2:ncol(xtab)) {
    xstats <- rbind(xstats,
      summary(as.freqtab(xtab[, c(1, i)])))
    ystats <- rbind(ystats,
      summary(as.freqtab(ytab[, c(1, i)])))
  }
  rownames(xstats) <- rownames(ystats) <-
    colnames(xtab)[-1]
  xvstats <- yvstats <- NULL
  if(object$design == "nonequivalent groups") {
    for(i in 2:ncol(xvtab)) {
      xvstats <- rbind(xvstats,
        summary(as.freqtab(xvtab[, c(1, i)])))
      yvstats <- rbind(yvstats,
        summary(as.freqtab(yvtab[, c(1, i)])))
    }
    rownames(xvstats) <- rownames(yvstats) <-
      colnames(xvtab)[-1]
  }
  
  out <- c(object[names(object) %in% c("name", "type", "method",
    "design", "ident", "ws", "smoothmethod")],
    freqtab = list(list(x = xtab, y = ytab,
      yx = yxtab)),
    stats = list(list(x = xstats, y = ystats,
      yx = yxstats)))
  if(!is.null(object$bootstraps))
    out$error <- summary(object$bootstraps)
  if(object$design == "nonequivalent groups") {
    out$freqtab <- c(out$freqtab, list(xv = xvtab,
      yv = yvtab))
    out$stats <- c(out$stats, list(xv = xvstats,
      yv = yvstats))
  }  
  class(out) <- "summary.equate"
  
  return(out)
}

# @describeIn equate Summary method for \dQuote{"equate.list"} objects.
#' @rdname equate
#' @export
summary.equate.list <- function(object, ...) {
  
  lapply(object, summary.equate)
}

#' @export
print.summary.equate <- function(x, ...) {
  
  cat("\n")
  cat(x$name, "\n\n")
  cat("Design:", x$design, "\n\n")
  if(x$type == "equipercentile") {
    cat("Smoothing Method: ")
    cat(switch(x$smoothmethod,
      bump = "adjusted frequency presmoothing",
      average = "average frequency presmoothing",
      loglinear = "loglinear presmoothing", "none"),
      "\n\n")
  }
  if(x$design == "nonequivalent groups" &&
      x$type != "composite" &&
      x$method != "chained")
    cat("Synthetic Weighting for x:",
      x$ws, "\n\n")
  lstats <- lapply(x$stats, function(y)
    t(round(y, 3)))
  out <- data.frame(matrix(unlist(lstats),
    ncol = ncol(x$stats[[1]]), byrow = TRUE))
  colnames(out) <- colnames(x$stats[[1]])
  scales <- names(unlist(lapply(lstats, colnames)))
  scales <- gsub("1|2|3|4", "", scales)
  dists <- unlist(lapply(lstats, colnames))
  dists <- gsub("(etic)|(erved)", "", dists)
  rownames(out) <- paste(scales, dists, sep = ".")
  cat("Summary Statistics:\n")
  print(out)
  cat("\n")
  
  if(!is.null(x$error)) {
    cat("Mean Error:\n")
    print(round(x$error, 4), row.names = FALSE)
    cat("\n")
  }
  
  invisible(x)
}
