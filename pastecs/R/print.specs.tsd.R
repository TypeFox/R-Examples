"print.specs.tsd" <-
function(x, ...) {
	print(unlist(unclass(x)), quote=FALSE)
	invisible(x)
}
