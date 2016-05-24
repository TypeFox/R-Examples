#' @method qread matrix
#' @export
qread.matrix <- function(
	file, type, header=TRUE, sep="\t", row.names=1, comment.char="",
	check.names=FALSE, dup.rm=FALSE, ...
) {
	if (dup.rm && !is.null(row.names)) {
		x <- read.table(file,
			header=header, sep=sep, comment.char=comment.char,
			check.names=check.names, ...);
		r <- x[,row.names];
		idx <- !duplicated(r);
		x <- as.matrix(x[idx, 2:ncol(x)]);
		rownames(x) <- r[idx];
		x
	} else {
		as.matrix(read.table(file,
			header=header, sep=sep, row.names=row.names, comment.char=comment.char,
			check.names=check.names, ...))
	}
}

#' @method qwrite matrix
#' @export
qwrite.matrix <- function(
	x, file, type, quote=FALSE, sep="\t", row.names=TRUE, col.names=NA, ...
) {
	if (!(is.matrix(x) || is.vector(x))) {
		stop("x must be a matrix or a vector");
	}
	write.table(x, file,
		quote=quote, sep=sep, row.names=row.names, col.names=col.names, ...)
}

#' @method qread mtx
#' @export
qread.mtx <- qread.matrix;

#' @method qwrite mtx
#' @export
qwrite.mtx <- qwrite.matrix;

#' @method qread dat
#' @export
qread.dat <- function(file, type, header=FALSE, sep="", comment.char="", ...) {
	as.matrix(read.table(file,
		header=header, sep=sep, comment.char=comment.char, ...));
}

#' @method qwrite dat
#' @export
qwrite.dat <- function(
	x, file, type, quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE, ...
) {
	if (!(is.matrix(x) || is.vector(x))) {
		stop("x must be a matrix or a vector");
	}
	write.table(x, file,
		quote=quote, sep=sep, row.names=row.names, col.names=col.names, ...)
}

