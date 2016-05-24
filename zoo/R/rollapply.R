
# data is a zoo object or a plain vector or matrix
#
# width is 
# - a list of integer vectors representing offsets or a plain vector of widths.
#   There is one per time point or its recycled if too short.  recycling uses
#   by= argument if length(width) is 1; otherwise, by is ignored.
#   If width represents widths then they are turned into offsets using align.
#
# If we are at 5th time of data and width[[5]] is c(-2,-1,0) then FUN is applied
#   to positions i + width[[i]] = 5 + c(-2,-1,0)  = 3:5 of z (so in terms of
#   a width specification it would be the same as width = 3, align = "right").
#
# Therefore we have the following transformations:
#   widths are converted to offsets which are converted to positions.  
#   The offsets are the components of width and the 
#   positions are i+width[[i]] after partial processing.  partial can be:
#   - logical.  FALSE means that all offets must exist or else no result is 
#     produced for that time point.  TRUE means that at least one offset must 
#     exist.
#   - numeric.  The minimum number of offsets that must exist.  If < 0 then
#     all elements of the offset must exist.  Note that TRUE corresponds to 1
#     and FALSE correspoinds to -1.  These are the two most common values.
#
# For points that are not computed they are filled in with fill.  fill has 
#   three elements and is recycled if too short.   fill = NULL is the default.
#   The elements represent what to fill the left points, interior points and
#   right points.  NULL causes no filling and "extend" causes the first or
#   last point to be repeated or interior points to be linearly approximated.

# wrapper around rollapply which defaults to align = "right"
rollapplyr <- function(..., align = "right") {
	rollapply(..., align = align)
}

rollapply <- function(data, ...) UseMethod("rollapply")

rollapply.default <- function(data, ...) {
	coredata(rollapply(zoo(data), ...))
}

rollapply.ts <- function(data, ...) {
	as.ts(rollapply(as.zoo(data), ...))
}

rollapply.zoo <- function(data, width, FUN, ..., by = 1, 
	by.column = TRUE, fill = if (na.pad) NA, na.pad = FALSE,
	partial = FALSE, align = c("center", "left", "right")) {

	if (!missing(na.pad)) {
		warning("na.pad argument is deprecated")
	}

    if (is.vector(width) && !is.list(width) && length(width) == 1 &&
		by.column && length(by) == 1 && by == 1 && (missing(partial) | identical(partial, FALSE)) &&
		length(list(...)) < 1 && length(sw <- deparse(substitute(FUN))) == 1) {
		  if (sw == "mean" && all(!is.na(data))) {
				return(rollmean(data, width, fill = fill, align = align))
		  } else if (sw == "median" && width %% 2 == 1) {
				return(rollmedian(data, width, fill = fill, align = align))
	      } else if (sw == "max") {
				return(rollmax(data, width, fill = fill, align = align))
	      }
	}
	FUN <- match.fun(FUN)

	if (by.column && length(dim(data)) == 2) {
		z <- do.call(merge,
			lapply(1:NCOL(data), function(j)
				rollapply(data[, j, drop = TRUE], width = width, FUN = FUN, ...,
					by = by, by.column = by.column, fill = fill,
					partial = partial, align = align)
			)
		)
		if (NCOL(data) == 1) dim(z) <- c(length(z), 1)
		colnames(z) <- if (NCOL(z) == NCOL(data)) colnames(data)
		return(z)
	}

	if (is.logical(partial)) partial <- if (partial) 1 else -1

	# convert widths to offsets using align
	align <- match.arg(align)

	if (!is.list(width)) width <- lapply(width, function(w) {
			if (align == "right") seq(to = 0, length.out = w)
			else if (align == "center") seq(to = floor(w/2), length.out = w)
			else seq(from = 0, length.out = w)
	})
	# recycle width (using by if length(width) == 1)
	width <- if (length(width) == 1) {
		w <- rep(list(NULL), NROW(data))
		start.at <- if (partial < 0) max(-min(width[[1]]), 0) + 1 else 1
		replace(w, seq(start.at, NROW(data), by = by), width)
	} else rep(width, length.out = NROW(data))

	f <- if (is.null(dim(data))) {
		# undimensioned
		#
		# if FUN is to be evaluated at offsets for the ith point then calculate
		# positions, do partial processing and apply FUN
		function(i, offsets, data, ...) { 
			if (is.null(offsets)) return(NULL)
			posns <- i + offsets
			ix <- posns >= 1 & posns <= NROW(data)
			if (partial < 0) {
				if (all(ix)) FUN(data[posns], ...)
			} else if (sum(ix) >= partial) {
				FUN(data[replace(posns, !ix, 0)], ...)
			}
		}
    } else {
		# dimensioned
		#
		# same f as in TRUE leg except data[.] becomes data[.,]
		function(i, offsets, data, ...) { 
			if (is.null(offsets)) return(NULL)
			posns <- i + offsets
			ix <- posns >= 1 & posns <= NROW(data)
			if (partial < 0) {
				if (all(ix)) FUN(data[posns,], ...)
			} else if (sum(ix) >= partial) {
				FUN(data[replace(posns, !ix, 0),], ...)
			}
		}
	}

	dat <- mapply(f, seq_along(time(data)), width, 
		MoreArgs = list(data = coredata(data), ...), SIMPLIFY = FALSE) 
		
	ix <- !sapply(dat, is.null) # integer indexes of non-nulls

	if (!missing(fill) || !missing(na.pad)) {

		# replace NULLs with NAs
		dat <- lapply(dat, function(x) if (is.null(x)) NA else x)

		# construct zoo object
		dat <-
		if (max(sapply(dat, length)) > 1)
			zoo(do.call("rbind", dat), index(data), attr(data, "frequency"))
		else
			zoo(do.call("c", dat), index(data), attr(data, "frequency"))

		# perform filling
		dat <- na.fill(dat, fill, ix)

	} else {

		# construct zoo object removing points corresponding to NULL
		dat <- if (max(sapply(dat, length)) > 1)
			zoo(do.call("rbind", dat), index(data)[ix], attr(data, "frequency"))
		else
			zoo(do.call("c", dat), index(data)[ix], attr(data, "frequency"))
	}

	dat
}
