"print.turnpoints" <-
function(x, ...) {
	cat("Turning points for:", x$data, "\n\n")
	cat("nbr observations  :", x$n, "\n")
	cat("nbr ex-aequos     :", sum(x$exaequos), "\n")
	if (x$firstispeak) {
		cat("nbr turning points:", x$nturns, "(first point is a peak)\n")
	} else {
		cat("nbr turning points:", x$nturns, "(first point is a pit)\n")
	}
	cat("E(p) =", 2 / 3 * (x$n - 2), "Var(p) =", (16 * x$n - 29) / 90, "(theoretical)\n")
	cat("\n")
	invisible(x)
}
