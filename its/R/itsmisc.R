
as.list.its <- function(x,...) {
	result <- if (ncol(x) == 1) list(x)
	else lapply(seq(length = ncol(x)), function(i) x[,i])
	names(result) <- colnames(x)
	result
}

as.data.frame.its <- function(x, row.names = NULL, optional = FALSE, ...) {
	result <- lapply(as.list(x), I)
	result$check.names <- result$check.rows <- optional
	do.call(data.frame, result)
}
