"print.turnogram" <-
function(x, ...) {
	cat(x$type, "turnogram for:", x$data, "\n\n")
	cat("options      :", x$fun, "/", x$proba, "\n")
	cat("intervals    :", min(x$interval), "..", max(x$interval), x$units.text, "/ step =", x$interval[2] - x$interval[1], "\n")
	cat("nbr of obs.  :", max(x$n), "..", min(x$n), "\n")
	maxinfo <- max(x$info)
	pos <- x$info == maxinfo
	cat("max. info.   : ", max(x$info), " at interval ", x$interval[pos], " (P = ", 2^-abs(max(x$info)), ": ", x$turns[pos], " turning points for ", x$n[pos], " observations)", "\n", sep="")
	cat("extract level: ", x$level, " (", x$units.text, ")\n\n", sep="")
	invisible(x)
}
