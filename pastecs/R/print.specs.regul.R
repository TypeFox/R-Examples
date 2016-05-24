"print.specs.regul" <-
function(x, ...) {
	print(unlist(unclass(x)), quote=FALSE)
	invisible(x)
}
