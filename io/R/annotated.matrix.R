
#' @method qread annotated.matrix
#' @export
qread.annotated.matrix <- function(
	file, type, annot.cols=1, header=TRUE, sep="\t", comment.char="", quote="", check.names=FALSE, stringsAsFactors=FALSE, ...
) {
	x <- read.table(file,
		header=header, sep=sep, comment.char=comment.char, quote=quote,
		check.names=check.names, stringsAsFactors=stringsAsFactors, ...);

	if (length(annot.cols) == 1) {
		# Expand annot.cols into a vector
		annot.cols <- 1:annot.cols;
	}

	# Extract meta information and data
	meta <- data.frame(x[, annot.cols], check.names=check.names);
	colnames(meta) <- colnames(x)[annot.cols];
	data <- as.matrix(x[, -annot.cols]);
	rownames(data) <- meta[, 1];
	
	structure(list(meta=meta, data=data), class="annotated.matrix")
}

#' @method qwrite annotated.matrix
#' @export
qwrite.annotated.matrix <- function(
	x, file, type, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE, ...
) {
	# Construct data.frame
	d <- data.frame(x$meta, x$data, check.names=FALSE);

	write.table(d, file,
		quote=quote, sep=sep, row.names=row.names, col.names=col.names, ...)
}

#' @method qread amtx
#' @export
qread.amtx <- function(file, type, annot.cols=1, ...) {
	qread.annotated.matrix(file, annot.cols=annot.cols)
}

#' @method qwrite amtx
#' @export
qwrite.amtx <- function(x, file, type, ...) {
	qwrite.annotated.matrix(x, file)
}
