"print.summary.turnogram" <-
function(x, ...) {
	cat(x$type, "turnogram for:", x$data, "\n\n")
	cat("options      :", x$fun, "/", x$proba, "\n")
	cat("call         :", deparse(x$call), "\n")
	maxinfo <- max(x$info)
	pos <- x$info == maxinfo
	cat("max. info.   : ", max(x$info), " at interval ", x$interval[pos], " (P = ", 2^-abs(max(x$info)), ": ", x$turns[pos], " turning points for ", x$n[pos], " observations)", "\n", sep="")
	cat("extract level: ", x$level, " (", x$units.text, ")\n\n", sep="")
	if (x$type == "Complete") {
		data <- list(interval=x$interval, n=x$n, turns=x$turns, turns.min=x$turns.min, turns.max=x$turns.max, info=x$info, info.min=x$info.min, info.max=x$info.max)
	} else {
		data <- list(interval=x$interval, n=x$n, turns=x$turns, info=x$info)
	}
	print(as.data.frame(data))
	invisible(x)
}
