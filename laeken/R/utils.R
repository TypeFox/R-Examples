# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

# TODO: error handling

#' Utility functions for indicators on social exclusion and poverty
#'
#' Test for class, print and take subsets of indicators on social exclusion and
#' poverty.
#'
#' @name utils
#'
#' @param x for \code{is.xyz}, any object to be tested.  The \code{print} and
#' \code{subset} methods are called by the generic functions if an object of the
#' respective class is supplied.
#' @param years an optional numeric vector giving the years to be extracted.
#' @param strata an optional vector giving the domains of the breakdown to be
#' extracted.
#' @param \dots additional arguments to be passed to and from methods.
#'
#' @return \code{is.indicator} returns \code{TRUE} if \code{x} inherits from
#' class \code{"indicator"} and \code{FALSE} otherwise.
#'
#' \code{is.arpr} returns \code{TRUE} if \code{x} inherits from class
#' \code{"arpr"} and \code{FALSE} otherwise.
#'
#' \code{is.qsr} returns \code{TRUE} if \code{x} inherits from class
#' \code{"qsr"} and \code{FALSE} otherwise.
#'
#' \code{is.rmpg} returns \code{TRUE} if \code{x} inherits from class
#' \code{"rmpg"} and \code{FALSE} otherwise.
#'
#' \code{is.gini} returns \code{TRUE} if \code{x} inherits from class
#' \code{"gini"} and \code{FALSE} otherwise.
#'
#' \code{is.gini} returns \code{TRUE} if \code{x} inherits from class
#' \code{"gini"} and \code{FALSE} otherwise.
#'
#' \code{print.indicator}, \code{print.arpr} and \code{print.rmpg} return
#' \code{x} invisibly.
#'
#' \code{subset.indicator}, \code{subset.arpr} and \code{subset.rmpg} return a
#' subset of \code{x} of the same class.
#'
#' @seealso \code{\link{arpr}}, \code{\link{qsr}}, \code{\link{rmpg}},
#' \code{\link{gini}}, \code{\link{gpg}}
#'
#' @references
#' A. Alfons and M. Templ (2013) Estimation of Social Exclusion Indicators
#' from Complex Surveys: The \R Package \pkg{laeken}.  \emph{Journal of
#' Statistical Software}, \bold{54}(15), 1--25.  URL
#' \url{http://www.jstatsoft.org/v54/i15/}
#'
#' @keywords survey
#'
#' @examples
#' data(eusilc)
#'
#' # at-risk-of-poverty rate
#' a <- arpr("eqIncome", weights = "rb050",
#'     breakdown = "db040", data = eusilc)
#' print(a)
#' is.arpr(a)
#' is.indicator(a)
#' subset(a, strata = c("Lower Austria", "Vienna"))
#'
#' # quintile share ratio
#' q <- qsr("eqIncome", weights = "rb050",
#'     breakdown = "db040", data = eusilc)
#' print(q)
#' is.qsr(q)
#' is.indicator(q)
#' subset(q, strata = c("Lower Austria", "Vienna"))
#'
#' # relative median at-risk-of-poverty gap
#' r <- rmpg("eqIncome", weights = "rb050",
#'     breakdown = "db040", data = eusilc)
#' print(r)
#' is.rmpg(r)
#' is.indicator(r)
#' subset(r, strata = c("Lower Austria", "Vienna"))
#'
#' # Gini coefficient
#' g <- gini("eqIncome", weights = "rb050",
#'     breakdown = "db040", data = eusilc)
#' print(g)
#' is.gini(g)
#' is.indicator(g)
#' subset(g, strata = c("Lower Austria", "Vienna"))
#'

NULL


## constructors
# class "indicator"
constructIndicator <- function(value, valueByStratum = NULL,
		varMethod = NULL, var = NULL, varByStratum = NULL, ci = NULL,
		ciByStratum = NULL, alpha = NULL, years = NULL, strata = NULL) {
	# construct and assign class
	x <- list(value=value, valueByStratum=valueByStratum, varMethod=varMethod,
			var=var, varByStratum=varByStratum, ci=ci, ciByStratum=ciByStratum,
			alpha=alpha, years=years, strata=strata)
	class(x) <- "indicator"
	# return object
	return(x)
}

# class "arpr"
constructArpr <- function(..., p = 0.6, threshold) {
	x <- constructIndicator(...)     # call constructor of superclass
	x$p <- p                         # set specific
	x$threshold <- threshold         # attributes
	class(x) <- c("arpr", class(x))  # assign class
	return(x)                        # return result
}

# class "qsr"
constructQsr <- function(...) {
	x <- constructIndicator(...)    # call constructor of superclass
	class(x) <- c("qsr", class(x))  # assign class
	return(x)                       # return result
}

# class "gpg"
constructGpg <- function(...) {
	x <- constructIndicator(...)    # call constructor of superclass
	class(x) <- c("gpg", class(x))  # assign class
	return(x)                       # return result
}

# class "rmrpg"
constructRmpg <- function(..., threshold) {
	x <- constructIndicator(...)     # call constructor of superclass
	x$threshold <- threshold         # set specific attributes
	class(x) <- c("rmpg", class(x))  # assign class
	return(x)                        # return result
}

# class "gini"
constructGini <- function(...) {
	x <- constructIndicator(...)    # call constructor of superclass
	class(x) <- c("gini", class(x))  # assign class
	return(x)                       # return result
}

# class "prop"
constructProp <- function(...) {
  x <- constructIndicator(...)    # call constructor of superclass
  class(x) <- c("prop", class(x))  # assign class
  return(x)                       # return result
}


## test for class

#' @rdname utils
#' @export
is.indicator <- function(x) inherits(x, "indicator")

#' @rdname utils
#' @export
is.arpr <- function(x) inherits(x, "arpr")

#' @rdname utils
#' @export
is.qsr <- function(x) inherits(x, "qsr")

#' @rdname utils
#' @export
is.rmpg <- function(x) inherits(x, "rmpg")

#' @rdname utils
#' @export
is.gini <- function(x) inherits(x, "gini")

#' @rdname utils
#' @export
is.prop <- function(x) inherits(x, "prop")

#' @rdname utils
#' @export
is.gpg <- function(x) inherits(x, "gpg")

## print
# class "indicator"
#' @rdname utils
#' @method print indicator
#' @export
print.indicator <- function(x, ...) {
	cat("Value:\n")
	print(x$value, ...)
	if(!is.null(x$var)) {
		cat("\nVariance:\n")
		print(x$var, ...)
	}
	if(!is.null(x$ci)) {
		cat("\nConfidence interval:\n")
		print(x$ci, ...)
	}
	if(!is.null(x$valueByStratum)) {
		cat("\nValue by domain:\n")
		print(x$valueByStratum, ...)
	}
	if(!is.null(x$varByStratum)) {
		cat("\nVariance by domain:\n")
		print(x$varByStratum, ...)
	}
	if(!is.null(x$varByStratum)) {
		cat("\nConfidence interval by domain:\n")
		print(x$ciByStratum, ...)
	}
	invisible(x)
}

# class "arpr"
#' @rdname utils
#' @method print arpr
#' @export
print.arpr <- function(x, ...) {
	print.indicator(x, ...)
	cat("\nThreshold:\n")
	print(x$threshold, ...)
	invisible(x)
}

# class "rmpg"
#' @rdname utils
#' @method print rmpg
#' @export
print.rmpg <- function(x, ...) {
	print.indicator(x, ...)
	cat("\nThreshold:\n")
	print(x$threshold, ...)
	invisible(x)
}

# class "minAMSE"
#' @export
print.minAMSE <- function(x, ...) {
	cat("Optimal k:\n")
	print(x$kopt, ...)
	cat("\nScale parameter:\n")
	print(x$x0, ...)
	cat("\nShape parameter:\n")
	print(x$theta, ...)
	invisible(x)
}


## subsets of indicators
# class "indicator"
#' @rdname utils
#' @method subset indicator
#' @export
subset.indicator <- function(x, years = NULL, strata = NULL, ...) {
	# initializations
	haveYears <- length(x$years) > 1
	haveVar <- !is.null(x$varMethod)
	haveStrata <- length(x$strata) > 1
	subsetYears <- haveYears && !is.null(years)
	subsetStrata <- haveStrata && !is.null(strata)
	# error handling
	if(subsetYears && !is.numeric(years)) {
		stop("'years' must be of type numeric")
	}
	if(subsetStrata && !is.character(strata)) {
		stop("'years' must be of type character")
	}
	# extract years from overall values (if available and requested)
	if(subsetYears) {
		ys <- as.character(years)
		x$value <- x$value[ys]
		if(haveVar) {
			x$var <- x$var[ys]
			x$ci <- x$ci[ys, , drop=FALSE]
		}
		x$years <- years  #set new years
	}
	# extract strata from overall values (if available and requested)
	if(subsetStrata || (haveStrata && subsetYears)) {
		n <- nrow(x$valueByStratum)
		if(subsetStrata) keepStrata <- x$valueByStratum$stratum %in% strata
		else keepStrata <- rep.int(TRUE, n)
		if(subsetYears) keepYears <- x$valueByStratum$year %in% years
		else keepYears <- rep.int(TRUE, n)
		keep <- keepStrata & keepYears
		x$valueByStratum <- x$valueByStratum[keep, , drop=FALSE]
		if(haveVar) {
			x$varByStratum <- x$varByStratum[keep, , drop=FALSE]
			x$ciByStratum <- x$ciByStratum[keep, , drop=FALSE]
		}
		x$strata <- strata  # set new strata
	}
	# return result
	return(x)
}

# class "arpr"
# TODO: allow for subsetting by threshold percentage
#' @rdname utils
#' @method subset arpr
#' @export
subset.arpr <- function(x, years = NULL, strata = NULL, ...) {
	haveYear <- length(x$years) > 1
	x <- subset.indicator(x, years, strata, ...)  # call method for superclass
  # subset threshold (if requested and available for multiple years)
	if(haveYear && !is.null(years)) {
		x$threshold <- x$threshold[as.character(years)]
	}
	# return result
	return(x)
}

# class "rmpg"
#' @rdname utils
#' @method subset rmpg
#' @export
subset.rmpg <- function(x, years = NULL, strata = NULL, ...) {
	haveYear <- length(x$years) > 1
	x <- subset.indicator(x, years, strata, ...)  # call method for superclass
	# subset threshold (if requested and available for multiple years)
	if(haveYear && !is.null(years)) {
		x$threshold <- x$threshold[as.character(years)]
	}
	# return result
	return(x)
}


## other utility functions

# get argument names of a function
argNames <- function(fun, removeDots = TRUE) {
	nam <- names(formals(fun))
	if(removeDots) nam <- setdiff(nam, "...")
	nam
}

# check percentages for the ARPT
checkP <-function(p) {
  if(is.numeric(p)) {
    keep <- !is.na(p) & p >= 0 & p <= 1
    p <- p[keep]
  } else p <- numeric()
  if(length(p) == 0) stop("'p' must contain numeric values in [0,1]")
  p
}

# get labels for percentages of the ARPT
getPLabels <-function(p) paste(signif(p*100), "%", sep="")
