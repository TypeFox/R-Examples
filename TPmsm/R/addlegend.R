addlegend <- function(legend.pos, curvlab, col=par("col"), lty=0, legend.bty, ...) {
	if ( is.list(legend.pos) ) legend.pos <- unlist(legend.pos)
	if (length(legend.pos) == 1) {
		xx <- legend.pos
		yy <- NULL
	}
	if (length(legend.pos) == 2) {
		xx <- legend.pos[1]
		yy <- legend.pos[2]
	}
	args <- list(...)
	ii <- pmatch( names(args), names(formals("legend")[-charmatch( "bty", names( formals("legend") ) )]) )
	do.call( "legend", c(list(xx, yy, curvlab, col=col, lty=lty, bty=legend.bty), args[!is.na(ii)]) )
	invisible()
}
