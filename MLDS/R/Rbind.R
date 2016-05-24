Rbind <- function(...) {
	.Deprecated("rbind", "MLDS")
	allargs <- list(...)
	dd <- do.call("rbind.data.frame", allargs)
	attr(dd, "invord") <- do.call("c", lapply(allargs,
		function(x) attr(x, "invord")))
#		attr(dd, "stimulus") <- attr(, "stimulus")
	class(dd) <- c("mlds.df", "data.frame")
	dd
	}