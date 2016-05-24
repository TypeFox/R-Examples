"print.summary.regul" <-
function(x, ...) {
	if (is.null(names(x$y))) {		# Only one vector is regulated
		cat("Regulation using method :", x$specs$methods, "\n")
		dat <- as.data.frame(cbind(x$x, x$y, x$match.dist))
		names(dat) <- c("Time", "Series", "Regul")
	} else {						# y is a data.frame
		cat("Regulation of, by \"method\" :\n")
		methods <- x$specs$methods
		names(methods) <- names(x$y)
		print(methods)
		dat <- as.data.frame(cbind(x$x, x$y, x$match.dist))
		Names <- NULL
		nc <- ncol(x$y) + 2
		Names[1] <- "Time"
		Names[2:(nc-1)] <- names(x$y)
		Names[nc] <- "Match.obs"
		names(dat) <- Names
	}
	cat("\nArguments for \"methods\" :\n")
	args <- NULL
	args[1] <- x$specs$tol.type
	args[2] <- x$specs$tol
	args[3] <- x$specs$rule
	args[4] <- x$specs$f
	args[5] <- x$specs$periodic
	args[6] <- x$specs$window
	args[7] <- x$specs$split
	names(args) <- c("tol.type", "tol", "rule", "f", "periodic", "window", "split")
	print(args)
	if (x$specs$rule == 1) {
		cat("\n", sum(x$match.dist == 1/0), "interpolated values on", length(x$match.dist), "(", sum(x$match.dist == -1/0), "NAs padded at ends )\n")
	} else {			# We allowed extrapolation
		cat("\n", sum(x$match.dist == 1/0), "interpolated and", sum(x$match.dist == -1/0), "extrapolated values on", length(x$match.dist), "\n")
	}		# Rem: 1/0 stands for Inf and -1/0 stands for -Inf
	cat("\nTime scale :\n")
	tsp <- NULL
	if (length(x$tspar$start) == 1) {
		tsp[1] <- x$tspar$start[1]
	} else {
		tsp[1] <- x$tspar$start[1] + (x$tspar$start[2] - 1) * x$tspar$frequency
	}
	tsp[2] <- x$tspar$deltat
	tsp[3] <- x$tspar$frequency
	names(tsp) <- c("start", "deltat", "frequency")
	print(tsp)
	cat("Time units :", x$units, "\n\n")
	cat("call :", deparse(x$call), "\n\n")
	print(dat)
	invisible(x)
}
