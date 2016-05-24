#' Frequency Distribution Tables
#' 
#' Functions for creating and manipulating frequency tables of class
#' \dQuote{\code{freqtab}}.
#' 
#' \code{freqtab} creates a frequency table from a vector or \code{data.frame}
#' of scores. When the \code{items} argument is included, scores are assumed to
#' be item responses, which are summed to create total scores. The scores are
#' tabulated and stored as an array, with dimensions for each variable. Note
#' that in previous versions of the \dQuote{\code{freqtab}} class the frequency
#' table was stored as a \code{data.frame}. This is no longer the case.
#' Instead, the table is stored as an array and converted to a
#' \code{data.frame} when printed or manipulated with the \code{head} and
#' \code{tail} methods.
#' 
#' \code{as.data.frame} converts an object of class \dQuote{\code{freqtab}} to
#' \dQuote{\code{data.frame}}. \code{droplevels} returns \code{x} with any
#' unused factor levels, or levels with zero counts, removed.
#' 
#' When \code{x} is an object of class \dQuote{\code{table}}, \code{freqtab}
#' simply modifies the attributes and converts to class
#' \dQuote{\code{freqtab}}. In this case, \code{x} must already be structured
#' similar to a \dQuote{\code{freqtab}} object, with the first dimension
#' containing counts for total scores, and remaining dimensions containing
#' counts for one or more anchor tests.
#' 
#' \code{as.freqtab} converts a 'flat' contingency table (see
#' \code{\link{ftable}}) to class \dQuote{\code{freqtab}} with the appropriate
#' attributes. A flat contingency table is the \code{data.frame} version of a
#' \dQuote{\code{freqtab}} object, where the first column contains the total
#' score scale, the last column contains counts, and the columns in between
#' contain different anchor test score combinations. \code{is.freqtab} tests
#' for class \dQuote{\code{freqtab}}.
#' 
#' \code{scales} extracts the measurement scales for the variables specified in
#' \code{margin}, with \code{margin = 1} referring to the total score scale,
#' and subsequent margins referring to anchor tests. \code{margin} is a wrapper
#' for \code{\link{margin.table}}, which itself is a simple wrapper for summing
#' over marginal counts, i.e., \code{apply(x, margin, sum)}. And \code{margins}
#' returns the number of dimensions, i.e., score variables, in a frequency
#' table.
#' 
#' \code{design} is used to set the dimnames of the frequency table, with
#' \code{total1} and \code{total2} used with single and counterbalanced groups,
#' and \code{total} and \code{anchor}(s) used otherwise. \code{design} also
#' sets the design attribute, which is used in \code{\link{equate}}.
#' 
#' The main difference between the \dQuote{\code{freqtab}} class and other
#' tabulation classes, like \dQuote{\code{table}} and \dQuote{\code{ftable}},
#' is that the \code{dimnames}, i.e., the score scales, are required to be
#' numeric. This facilitates plotting with \code{\link{plot.freqtab}}, equating
#' with the \code{\link{equate}} function, and descriptive statistics with the
#' \code{\link{summary.freqtab}} and other methods.
#' 
#' @param x either an object (vector or \code{data.frame}) containing total
#' scores or item responses with which total scores will be calculated, or an
#' object inheriting from class \dQuote{\code{freqtab}}. In the \code{freqtab}
#' function, the first column in \code{x} must be the total test; any remaining
#' columns may contain anchor scores. See below for details.
#' @param scales list of vectors containing the score scales for each score scale
#' in \code{x}.
#' @param items list of vectors of column indices (numbers or column names)
#' with which total scores will be computed for \code{x} when it contains item
#' responses.
#' @param design the equating design used in data collection. For univariate
#' \code{x}, \code{design = "eg"} is assumed for equivalent groups. For
#' multivariate \code{x}, \code{design = "ng"} is assumed for nonequivalent
#' groups. Single-groups and counterbalanced designs must be specified with
#' \code{design = "sg"} and \code{design = "cb"}.
#' @param na.rm logical with default \code{TRUE} specifying whether or not
#' missing item responses should be ignored, i.e., treated as 0, when
#' calculating total scores.
#' @param row.names,optional arguments passed to \code{as.data.frame},
#' currently ignored.
#' @param drop logical, with default \code{FALSE}, indicating whether or not
#' unused factor levels, or score values with zero counts, should be dropped.
#' See below for details.
#' @param margin integer vector specifying the margin(s) over which frequencies
#' should be summed.
#' @param ... further arguments passed to or from other functions.
#' @return A table array with dimensions equal to the number of score scales.
#' In most cases, this will be a univariate or bivariate distribution, but
#' multivariate distributions are supported. \code{scales} and \code{margins}
#' return numeric vectors.
#' @author Anthony Albano \email{tony.d.albano@@gmail.com}
#' @seealso \code{\link{table}}, \code{\link{ftable}},
#' \code{\link{summary.freqtab}}, \code{\link{plot.freqtab}}
#' @keywords methods
#' @examples
#' 
#' # Univariate distribution with score scale
#' set.seed(2005)
#' x <- round(rnorm(1000, 100, 10))
#' head(freqtab(x, scales = 70:130))
#' 
#' # Existing frequency table converted to class "freqtab"
#' # The first score of zero, with zero counts, is dropped
#' head(as.freqtab(ACTmath[, 1:2], drop = TRUE))
#' 
#' # Bivariate distribution
#' # Reduce y to the anchor test margin (2)
#' ny <- freqtab(x = KBneat$y, scales = list(0:36, 0:12))
#' margin(ny, margin = 2)
#' 
#' # Summing scored item responses with the PISA data
#' attach(PISA)
#' r6items <- paste(items$itemid[items$clusterid == "r6"])
#' r7items <- paste(items$itemid[items$clusterid == "r7"])
#' pisa67 <- freqtab(students[students$book == 6, ],
#' 	items = list(r6items, r7items),
#' 	scales = list(0:16, 0:14))
#' detach(PISA)
#' 
#' # Scales for both margins
#' # Zero total score is unobserved
#' scales(pisa67, 1:2)
#' scales(droplevels(pisa67), 1:2)
#' 
#' @export freqtab
freqtab <- function(x, ...) UseMethod("freqtab")

#' @describeIn freqtab Default method for a \code{data.frame} of item responses,
#' a \code{data.frame} of total and anchor scores, or a vector of total scores.
#' @export
freqtab.default <- function(x, scales, items, design,
	na.rm = TRUE, ...) {

	x <- data.frame(x)
	if (!missing(items)) {
		if (!is.list(items))
			items <- list(items)
		x <- sapply(items, function(i) {
			if (is.factor(i))
				i <- as.character(i)
			rowSums(cbind(x[, i]), na.rm = na.rm)
			})
	}

	nx <- ncol(x)
	x <- as.data.frame(x,
		row.names = seq_len(nrow(x)))
	if (missing(scales))
		scales <- apply(x, 2, function(y)
			min(y):max(y))
	if (!is.list(scales))
		scales <- list(scales)
	for (i in 1:nx)
		x[, i] <- factor(as.character(x[, i]),
			levels = scales[[i]])
	
	out <- freqtab.table(table(x), design)

	return(out)
}

# Note drop will drop levels with zero counts across all
# other factors - you can't drop if x is a vector
# adding levels with zero counts is achieved by
# including those values in scales and setting drop
# to FALSE
# You can't drop current zeros and then add new ones

#' @rdname freqtab
#' @export
as.freqtab <- function(x, scales, design, drop = FALSE, ...) {
	
	x <- as.data.frame(x)
	nx <- ncol(x) - 1

	# If x is a vector of counts, scales must be supplied
	if (nx == 0) {
		dm <- sapply(scales, length)
		if (nrow(x) != prod(dm))
			stop("'scales' do not match dimensions of 'x'")
		out <- array(unlist(x), dm, dimnames = scales)
	}
	else {
		if (drop)
			x <- x[x[, nx + 1] > 0, ]
		# Set scales
		if (missing(scales)) {
			if (nx == 1)
				scales <- list(unique(x[, 1]))
			else if (nx > 1)
				scales <- lapply(x[, -(nx + 1)], function(y)
					as.character(sort(unique(y))))
			else
				stop("'scales' cannot be inferred from 'x'")
		}
		else if (!is.list(scales))
			scales <- list(scales)
		# Refactor like in xtabs
		xf <- lapply(1:nx, function(i) {
			factor(as.character(x[, i]),
				levels = scales[[i]])[, drop = drop]
		})
    	out <- tapply(x[, nx + 1], xf, sum)
	    out[is.na(out)] <- 0
	}

	return(freqtab.table(as.table(out), design))
}

is.freqtab <- function(x) {
	
	return(class(x)[1] == "freqtab")
}

#' @describeIn freqtab Method for \code{table}s.
#' @export
freqtab.table <- function(x, design, ...) {
	
	out <- c(x)
	attributes(out) <- attributes(x)
	nx <- margins(x)
	if (missing(design)) {
		if(nx == 1) design <- "eg"
		else if(nx > 1) design <- "ng"
	}
	design <- match.arg(design, c("sg", "cb", "eg", "ng"))
	if (design %in% c("sg", "cb")) {
		dnn <- c("total1", "total2")
	}
	else {
		dnn <- "total"
		if (nx == 2)
			dnn[2] <- "anchor"
		else if (nx > 2)
			dnn[2:nx] <- paste("anchor",
				1:(nx - 1), sep = "")
	}
	names(dimnames(out)) <- dnn
	attr(out, "design") <- design
	class(out) <- c("freqtab", "table")
	
	return(out)
}

design <- function(x) {
	if(is.freqtab(x)) out <- attributes(x)$design
  else out <- NA
  return(out)
}

is.design <- function(x, type = c("sg", "cb", "eg", "ng")) {
  type <- match.arg(type)
  return(ifelse(design(x) == type, TRUE, FALSE))
}

# @describeIn freqtab Convert \dQuote{\code{freqtab}} object to
# \dQuote{\code{data.frame}}.
#' @rdname freqtab
#' @export
as.data.frame.freqtab <- function(x, row.names = NULL,
	optional = FALSE, drop = FALSE, ...) {
	
	out <- as.data.frame.table(x, row.names = NULL,
		responseName = "count", stringsAsFactors = TRUE)
	out <- sapply(out, function(y)
		as.numeric(as.character(y)))
	if (drop)
		out <- out[out[, ncol(out)] > 0, ]

	return(as.data.frame(out))
}

# @describeIn freqtab \code{head} method for \dQuote{\code{freqtab}} objects.
#' @rdname freqtab
#' @importFrom utils head
#' @export
head.freqtab <- function(x, ...) {
	
	utils::head(as.data.frame(x), ...)
}

# @describeIn freqtab \code{tail} method for \dQuote{\code{freqtab}}.
#' @rdname freqtab
#' @importFrom utils tail
#' @export
tail.freqtab <- function(x, ...) {
	
	utils::tail(as.data.frame(x), ...)
}

#' @export
print.freqtab <- function(x, ...) {
	
	print(as.data.frame(x))
}

# @describeIn freqtab Extract the scale(s) from a \dQuote{\code{freqtab}} object.
#' @rdname freqtab
#' @export
scales <- function(x, margin = 1) {
	
	if (length(margin) == 1)
		return(as.numeric(dimnames(x)[[margin]]))
	else if (length(margin) > 1)
		return(lapply(dimnames(x)[margin], as.numeric))
}

# @describeIn freqtab Extract marginal distributions from a
# \dQuote{\code{freqtab}} object.
#' @rdname freqtab
#' @export
margin <- function(x, margin = 1) {
	
	if (any(!margin %in% seq(margins(x))))
		stop("misspecified margins")

	out <- margin.table(x, margin)
  attr(out, "design") <- design(x)
  return(out)
}

# @describeIn freqtab Find the number of margins in a
# \dQuote{\code{freqtab}} object.
#' @rdname freqtab
#' @export
margins <- function(x) {
	
	return(length(dim(x)))
}

# @describeIn freqtab Drop unused or unobserved levels in a
# \dQuote{\code{freqtab}} object.
#' @rdname freqtab
#' @export
droplevels.freqtab <- function(x, ...) {
	
	xd <- as.data.frame(x)[x > 0, ]
	return(as.freqtab(droplevels(xd)))
}
