`importance` <-
function (x) UseMethod("importance")

`importance.averaging` <-
function(x) return(x$importance)

`importance.model.selection` <-
function(x) {
	if(nrow(x) <= 1L) stop("argument consists of only one model")

	tt <- attr(x, "terms")
	z <- x[, tt, drop = FALSE]
	z <- !is.na(z[, !apply(apply(z, 2L, is.na), 2, all) & !(tt %in% attr(tt, "interceptLabel")),
		drop = FALSE])
	wt <- x[, "weight"]
	res <- apply(z, 2L, function(y) sum(wt[y]))
	o <- order(res, decreasing = TRUE)
	res <- res[o]
	attr(res, "n.models") <- colSums(z)[o]
	class(res) <- c("importance", "numeric") 
	return(res)
}


`print.importance` <-
function(x, ...) {
	print.default(format(matrix(c(
		format(ifelse(x < 0.01, "<0.01", zapsmall(x, 2L)), scientific = FALSE,
		justify = "r"), format(attr(x, "n.models"))), nrow = 2L, byrow = TRUE,
		dimnames = list(c("Importance:", "N containing models:"), names(x))),
		justify = "r"), quote = FALSE)
	invisible(x)
}

#function(x) return(apply(x[, attr(x, "terms")], 2L,
#	function(z) sum(x[, "weight"][!is.na(z)])))

`importance.default` <- function(x)
	model.avg(x)$importance