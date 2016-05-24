#' @method droplevels lfactor
#' @export
droplevels.lfactor <- function(x, ...) {
	labels <- attributes(x)$levels
	levels <- attributes(x)$llevels
	cookie_cutter <- labels %in% x
	labels <- labels[cookie_cutter]
	levels <- levels[cookie_cutter]
	return(lfactor(switchllevels(x), levels=levels, labels=labels))
	#attr(x, "levels") <- labels
	#attr(x, "llevels") <- levels
	#return(x)
}
