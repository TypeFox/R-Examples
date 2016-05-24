print.klcv <- function (x, digits = max(3, getOption("digits") - 3), ...){
	rho <- x$rho
	out_gdf <- x$gdf
	out_ll <- x$loglik
	out_klcv <- x$klcv
	min.klcv <- x$min.klcv
	rho.opt <- x$rho.opt
	tbl <- data.frame(rho, out_gdf, out_ll, out_klcv)
	names(tbl) <- c("rho", "gdf", "log-lik", "klcv")
	tbl.format <- format(tbl, digits = digits)
	print(tbl.format, print.gap = 2, quote = FALSE, row.names=FALSE, ...)
	cat("\nminimum KLCV = ", min.klcv, "at rho =", rho.opt, "(scale factor =", x$scale, ")\n\n")
	invisible(tbl)
}
