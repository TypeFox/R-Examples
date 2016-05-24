#' @method qread gmt
#' @export
qread.gmt <- function(file, type, ...) {
	sep <- "\t";

	x <- lapply(
		strsplit(readLines(file, warn=FALSE), sep),
		function(z) {
			list(meta=z[1:2], data=z[3:length(z)])
		}
	)

	# organize s.t. meta and data fields are at the top level	
	meta <- lapply(x, function(z) z$meta[2]);
	names(meta) <- unlist(lapply(x, function(z) z$meta[1]));
	data <- lapply(x, function(z) z$data);
	names(data) <- names(meta);

	structure(list(meta=meta, data=data), class="gene.sets")
}

#' @method qwrite gmt
#' @export
qwrite.gmt <- function(x, file, type, ...) {
	sep <- "\t";

	if (!inherits(x, "gene.sets")) {
		stop("x must be a gene.sets");
	}

	# serialize x into a character vector
	sets <- lapply(x$data, function(z) paste(z, collapse=sep));
	serialized <- paste(names(x$meta), unlist(x$meta), sets, sep=sep);

	writeLines(serialized, file)
}
