
dyn <- function(x) {
	class(x) <- unique(c("dyn", oldClass(x)))
	x
}
dyn <- dyn(dyn)

"$.dyn" <- function(x, fun) {
	e <- parent.frame()
	if (!identical(x, dyn)) return(eval(substitute(unclass(x)$fun),e))
	# fun <- as.name(fun)
	f <- substitute(function(...) {
		cl <- match.call()
		cl[[1]] <- as.name(fun)
		cl[[2]] <- as.call(list(as.name("dyn"), cl[[2]]))
		result <- eval.parent(cl)
		dyn(result)
	})
	eval(f, e)
}

