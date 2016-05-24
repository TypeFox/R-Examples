
# fill is a 3 component list or is coerced to one representing
# fill char to left of leftmost non-NA, fill character to interior NAs
#  and fill char to right of rightmost non-NA
# If component is "extend" then left or rightmost NA is extended or interior
#  NA is linearly interpolated
# If component is NULL then the corresponding NA is dropped.

na.fill <- function(object, fill, ...) UseMethod("na.fill")

na.fill.zoo <- function(object, fill, ix, ...) {

	if (length(dim(object)) == 2 && NCOL(object) > 1) {
		ixmiss <- missing(ix)
		L <- lapply(1:NCOL(object), 
				function(i) {
					if (ixmiss) ix <- !is.na(object[,i])
					na.fill(object[,i], fill, ix, ...)
				})
		out <- do.call("merge", c(L, all = FALSE))
		colnames(out) <- colnames(object)
		return(out)
	}

	if (missing(ix)) ix <- !is.na(object)

	if ((is.logical(ix) && any(ix)) || (!is.logical(ix) && length(ix))) {

		n <- length(object)
		# integer indexes for output points which are present
		wix <- if (is.logical(ix)) which(ix) else ix
		# min and max integer index
		wx.min <- head(wix, 1) 
		wx.max <- tail(wix, 1)
		# similar to wrng <- wx.min:wx.max
		wrng <- seq(wx.min, length.out = wx.max - wx.min + 1)

		# recycle to length 3
		fill <- rep(as.list(fill), length.out = 3)
		# we will be coercing fill values to the class of coredata(data).
		# This allows fill=c("extend", NA) to work even though NA is coerced to
		#  a character NA.
		as.cls <- if (is.numeric(coredata(object))) as.numeric else as.logical
		fill <- lapply(fill, function(x) if (is.character(x) &&
			pmatch(x, "extend", nomatch = 0)) "extend" else as.cls(x))
		# fill points on left
		if (length(fill[[1]]) > 0) 
			if (!is.null(fill[[1]])) object[seq_len(wx.min - 1)] <- 
				if (is.character(fill[[1]]) && fill[[1]] == "extend")
						object[[wx.min]] else fill[[1]]
		# fill intermediate points
		# - this is for zoo method, for zooreg method it would be possible to
		#   perform linear interpolation in proportion to time rather than
		#   in proportion to the integer index
		if (length(fill[[2]]) > 0) {
			if (is.character(fill[[2]]) && fill[[2]] == "extend") object[wrng] <- 
					# as.list(approx(wix, unlist(object[wix]), xout = wrng)$y)
					approx(wix, unlist(object[wix]), xout = wrng)$y
			else object[intersect(which(!ix), wrng)] <- fill[[2]]
		}
		# fill points on right
		if (length(fill[[3]]) > 0) 
			object[seq(wx.max + 1, length.out = n - wx.max)] <- 
				if (is.character(fill[[3]]) && fill[[3]] == "extend")
						object[[wx.max]] else fill[[3]]

		keep <- seq_len(n)
		if (length(fill[[1]]) == 0) keep <- unique(pmax(wx.min, keep))
		if (length(fill[[2]]) == 0) {
			wrng <- seq(wx.min, length.out = wx.max - wx.min + 1)
			keep <- setdiff(keep, intersect(which(!ix), wrng))
		}
		if (length(fill[[3]]) == 0) keep <- unique(pmin(wx.max, keep)) 
		return(object[keep, , drop = is.null(dim(object))])
	} else if(length(fill)) {
	  object[is.na(object)] <- fill
	  return(object)
	}
}

na.fill.default <- function(object, fill, ix, ...) {
	coredata(na.fill(zoo(object), fill, ix, ...))
}
	
na.fill.ts <- function(object, fill, ix, ...) {
	as.ts(na.fill(as.zoo(object), fill, ix, ...))
}

