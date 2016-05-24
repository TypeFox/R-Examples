	mplot <- function(x, y, id, data=parent.frame(), subset=NULL, add=FALSE, ...) {
#	plots y ~ x by id with data
#	x and y can be name or character
#	subset defines a subset of rows
#	add TRUE suppresses plot axes
#	... parameters where bg, cex, col, lty, lwd, pch can depend on id

#	save x y id
	mcall <- match.call()[-1]
	df <- as.data.frame(lapply(as.list(mcall[1:3]), function(z) {
		if (is.character(z)) with(data,	get(z, inherits=TRUE))
		else eval(z, envir = data, enclos = parent.frame())
	}))
	if (length(deparse(mcall)) == 1) {
    names(df) <- lapply(as.list(mcall[1:3]), function(z) {
  		if (is.character(z)) z else deparse(z)
  	})
	}

#	extract and save vector par args: bg cex col lty lwd pch
	if (length(dots <- match.call(expand.dots=FALSE)$...) > 0) {
		ARG <- lapply(as.list(dots), eval, envir = data, enclos = parent.frame())
		cnames <- names(ARG)[lapply(ARG, length) == nrow(df)]
		df[, cnames] <- ARG[cnames]
		ARG[cnames] <- NULL
	} else {
		ARG <- list()
		cnames <- NULL
	}

#	subset data
	subset <- eval(substitute(subset), data, parent.frame())
	if (!is.null(subset)) {
		if (length(subset) != nrow(df)) stop('subset wrong length for data')
		subset <- ifelse(is.na(df[, 1]) | is.na(df[, 2]), FALSE, subset)
		df <- df[subset, ]
	}
	if (nrow(df) == 0) stop("no data to plot")

#	plot axes if new graph
	if (!add) {
		if (!"xlab" %in% names(ARG)) ARG <- c(ARG, list(xlab=quote(names(df)[1])))
		if (!"ylab" %in% names(ARG)) ARG <- c(ARG, list(ylab=quote(names(df)[2])))
		type <- match(names(ARG), "type", 0)
		do.call("plot", c(list(x=df[, 1], y=df[, 2], type='n'), ARG[!type]))
	}

#	draw growth curves
	tt <- by(df, df[, 3], function(z) {
#	sort by x
		ox <- order(z[, 1])
#	restore vector ... args
		if (length(cnames) > 0) ARG[cnames] <- as.list(as.data.frame(z[ox, cnames], stringsAsFactors=FALSE))
#	lines(x, y, ...)
		do.call("lines", c(list(x=z[ox, 1], y=z[ox, 2]), ARG))
	})
}
