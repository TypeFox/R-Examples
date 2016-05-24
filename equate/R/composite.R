#' Composite Linking and Equating
#' 
#' This function creates a composite linking or equating as a combination of
#' two or more other linking or equating functions.
#' 
#' Composite linking and equating create a single linking or equating function
#' as a weighted combination of two or more other linking or equating
#' functions. See Holland and Strawderman (2011) for details.
#' 
#' @param x for the default method, \code{x} is a matrix of equating functions,
#' with one function per column. Otherwise, \code{x} is a list of equatings,
#' where each element is an object of class \dQuote{\code{\link{equate}}}.
#' @param wc vector of weights for creating the composite. \code{length(wc)}
#' should match either \code{ncol(x)} for the default method or
#' \code{length(x)}.
#' @param name an optional name, used to label the output. If missing, a name
#' will be created using \code{x}.
#' @param p integer specifying the type of circle used to define symmetry.
#' @param symmetric logical, with default \code{FALSE}, indicating whether or
#' not weights \code{wc} should be modified to create symmetric weights. Only
#' supported for composites of linear functions.
#' @param verbose logical, with default \code{TRUE}, indicating whether or not
#' full output should be returned. When \code{FALSE}, only the equated scores
#' are returned.
#' @param \dots further arguments passed to or from other functions.
#' @return For the default method, and when \code{verbose = FALSE}, a vector of
#' composite equated scores is returned. Otherwise, a list of equating output
#' is returned, including output for the composite and each function being
#' combined.
#' @author Anthony Albano \email{tony.d.albano@@gmail.com}
#' @seealso \code{\link{equate}}
#' @references Holland, P. W., and Strawderman, W. E. (2011). How to average
#' equating functions, if you must. In A. A. von Davier (Ed.), Statistical
#' models for test equating, scaling, and linking (pp. 89-107). New York, NY:
#' Springer.
#' @keywords methods
#' @examples
#' # See vignette("equatevignette") for additional examples
#' 
#' # Example from the equate help file, without the bootstrapping
#' # Random groups equating for (1) identity, (2) mean, 
#' # (3) linear, (4) equipercentile with loglinear
#' # smoothing, and (5) a composite of mean and identity
#' rx <- as.freqtab(ACTmath[, 1:2])
#' ry <- as.freqtab(ACTmath[, c(1, 3)])
#' set.seed(2007)
#' 
#' req1 <- equate(rx, ry, type = "i")
#' req2 <- equate(rx, ry, type = "m")
#' req3 <- equate(rx, ry, type = "l")
#' req4 <- equate(rx, ry, type = "e", smooth = "loglin",
#'   degrees = 3)
#' req5 <- composite(list(req1, req2), wc = .5, symmetric = TRUE)
#' 
#' # Compare equating functions
#' plot(req1, req2, req3, req4, req5[[1]], addident = FALSE)
#' 
#' @export composite
composite <- function(x, ...) UseMethod("composite")

#' @describeIn composite Default method for a matrix of equating functions,
#' one per column.
#' @export
composite.default <- function(x, wc, ...) return(x %*% wc)

#' @describeIn composite Create composite across functions in
#' \dQuote{\code{equate.list}} object.
#' @export
composite.equate.list <- function(x, wc, name, symmetric = FALSE,
	p = 1, verbose = TRUE, ...) {
	
	if(missing(wc))
		wc <- rep(1/length(x), length(x))
	if(symmetric) {
		if(!all(sapply(x, function(z) z$type) %in%
			c("identity", "mean", "linear", "general"))) {
				warning("all components must be linear to create ",
					"symmetric weights")
				wcs <- wc
				symmetric <- FALSE
			}
		else {
			slopes <- sapply(x, function(z) z$coef[2])
			wcs <- (wc*(1 + slopes^p)^-(1/p))/
					sum(wc*(1 + slopes^p)^-(1/p))
		}
	}
	else wcs <- wc
		
	yx <- composite.default(sapply(x, function(z) z$conc$yx),
		wcs)

	if(verbose) {
		if(missing(name))
			name <- paste("Composite Equating:",
				strsplit(x[[1]]$name, ": ")[[1]][2])
		out <- list(name = name, type = "composite",
			design = x[[1]]$design, x = x[[1]]$x, y = x[[1]]$y,
			concordance = data.frame(scale = x[[1]]$conc$scale,
				yx = yx), wc = wc, symmetric = symmetric)
		if(symmetric) out$wcs <- wcs
		out <- as.composite(out)
		out <- as.equate.list(c(list(out), x))
	}
	else out <- yx
	
	return(out)
}

#' @describeIn composite Create composite across functions in
#' \dQuote{\code{list}} object.
#' @export
composite.list <- function(x, wc, name,	symmetric = FALSE,
	p = 1, verbose = TRUE, ...) {
  
	if(!all(sapply(x, function(z) is.equate(z))))
		stop("all elements of 'x' must be class 'equate'")
	
	return(composite(as.equate.list(x), wc, name,
		symmetric, p, verbose, ...))
}

as.composite <- function(x) {
	
	class(x) <- c("composite", "equate")
	return(x)
}

is.composite <- function(x) {
	
	return(class(x)[1] == "composite")
}

#' @export
print.composite <- function(x, ...) {

	cat("\n")
	cat(x$name, "\n\n")
	cat("Design:", x$design, "\n\n")

	stats <- rbind(x = summary(margin(x$x)),
		y = summary(margin(x$y)),
		yx = summary(as.freqtab(cbind(x$concordance[, 2],
			c(margin(x$x))))))
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
