"print.summary.turnpoints" <-
function(x, ...) {
	cat("Turning points for:", x$data, "\n\n")
	cat("nbr observations  :", x$n, "\n")
	cat("nbr ex-aequos     :", sum(x$exaequos), "\n")
	if (x$firstispeak) {
		cat("nbr turning points:", x$nturns, "(first point is a peak)\n")
		typep <- c("peak", "pit")
	} else {
		cat("nbr turning points:", x$nturns, "(first point is a pit)\n")
		typep <- c("pit", "peak")
	}
	cat("E(p) =", 2 / 3 * (x$n - 2), "Var(p) =", (16 * x$n - 29) / 90, "(theoretical)\n")
	cat("\n")
	# construct the table summarizing all turning points
	typepts <- rep(typep, length.out=x$nturns)
	tablepts <- as.data.frame(list(point=x$tppos, type=typepts, proba=x$proba, info=x$info))
	print(tablepts)
	invisible(x)
}
