na.approx <- function(object, ...) UseMethod("na.approx")

na.approx.zoo <- function(object, x = index(object), xout, ..., na.rm = TRUE, along) {

    if (!missing(along)) {
        warning("along to be deprecated - use x instead")
        if (missing(x)) x <- along
    }

    missing.xout <- missing(xout) || is.null(xout)
    if (is.function(x)) x <- x(index(object))
    if (!missing.xout && is.function(xout)) xout <- xout(index(object))
    order.by <- if (missing.xout) index(object) else xout
    xout <- if (missing.xout) x else xout

    if (missing.xout || identical(xout, index(object))) {
        result <- object
    } else {
        object.x <- object
        if (!identical(class(x), class(xout))) {
            index(object.x) <- as.numeric(x)
            xout <- as.numeric(xout)
        } else {
            index(object.x) <- x
        }
        objectm <- merge(object.x, zoo(, xout))
        if (length(dim(objectm)) == 2) colnames(objectm) <- colnames(object)
        result <- window(objectm, index = xout)
    }
    result[] <- na.approx.default(object, x = x, xout = xout, na.rm = FALSE, ...)
    if ((!missing(order.by) && !is.null(order.by)) || !missing.xout) {
        index(result) <- order.by
    }

    if (na.rm) {
        result <- na.trim(result, is.na = "all")
    }

    result

}

na.approx.zooreg <- function(object, ...) {
    object. <- structure(object, class = setdiff(class(object), "zooreg"))
    as.zooreg(na.approx(object., ...))
}


na.approx.default <- function(object, x = index(object), xout = x, ..., na.rm = TRUE, maxgap = Inf, along) {

    if (!missing(along)) {
        warning("along to be deprecated - use x instead")
        if (missing(x)) x <- along
    }

    na.approx.vec <- function(x, y, xout = x, ...) {
        na <- is.na(y)
	if(sum(!na) < 2L) {
	    ## approx() cannot be applied here, hence simply:
	    yf <- rep.int(NA, length(xout))
	    if(any(!na)) {
	        if(x[!na] %in% xout) {
		    yf[xout == x[!na]] <- y[!na]
		}
	    }
	    return(yf)
	}
        yf <- approx(x[!na], y[!na], xout, ...)$y
        if (maxgap < length(y)) {
            ## construct a series like y but with only gaps > maxgap
            ## (actual values don't matter as we only use is.na(ygap) below)
            ygap <- .fill_short_gaps(y, seq_along(y), maxgap = maxgap)
            ## construct y values at 'xout', keeping NAs from ygap
            ## (using indexing, as approx() does not allow NAs to be propagated)
            ix <- approx(x, seq_along(y), xout, ...)$y
            yx <- ifelse(is.na(ygap[floor(ix)] + ygap[ceiling(ix)]), NA, yf)
            yx
        } else {
            yf
        }
    }

    if (!identical(length(x), length(index(object)))) {
        stop("x and index must have the same length")
    }
    x. <- as.numeric(x)
    if (missing(xout) || is.null(xout)) xout <- x.
    xout. <- as.numeric(xout)
    object. <- coredata(object)

    result <- if (length(dim(object.)) < 2) {
        na.approx.vec(x., coredata(object.), xout = xout., ...)
    } else {
        apply(coredata(object.), 2, na.approx.vec, x = x., xout = xout., ...)
    }

    if (na.rm) {
        result <- na.trim(result, is.na = "all")
    }

    result

}

na.approx.ts <- function(object, ...) {
    as.ts(na.approx(as.zoo(object), ...))
}

## x = series with gaps
## fill = same series with filled gaps
.fill_short_gaps <- function(x, fill, maxgap) {
    if (maxgap <= 0)
        return(x)
    if (maxgap >= length(x))
        return(fill)
    naruns <- rle(is.na(x))
    naruns$values[naruns$lengths > maxgap] <- FALSE
    naok <- inverse.rle(naruns)
    ifelse(naok, fill, x)
}
