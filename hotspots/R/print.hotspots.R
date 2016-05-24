print.hotspots <-
function(x, ...) {
ty <- "hot spot"
if (!is.null(x$u)) ty <- "outlier"
if (!is.null(x$positive.cut)) {
	cat(paste(ty, "cutoff (positive): \n"))
	print(x$positive.cut) }
if (!is.null(x$negative.cut)) {
	cat(paste(ty, "cutoff (negative): \n"))
	print(x$negative.cut) }}

