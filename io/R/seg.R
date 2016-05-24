#' @method qread seg
#' @export
qread.seg <- function(file, type, ...) {

	# open a file connection
	if (is.character(file)) {
		f <- base::file(file, "rt");
	} else {
		f <- file;
	}

	header <- scan(f, character(), nlines=1, comment.char="#", quiet=TRUE);

	# continue reading the file connection
	# assume first four columns are:
	# name (character), chromosome (character), start (integer), end (integer)
	x <- read.table(f, sep="\t", header=FALSE,
		colClasses=c("character", "character", "integer", "integer",
			rep("numeric", length(header)-4)), ...);
	colnames(x) <- header;

	if (is.character(file)) {
		close(f);
	}

	x
}

#' @method qwrite seg
#' @export
qwrite.seg <- function(x, file, type, ...) {
	if (!is.data.frame(x)) {
		stop("x must be a data.frame");
	}
	write.table(x, file, sep="\t", row.names=FALSE, col.names=TRUE,
		quote=FALSE, na="")
}
