#' @method qread gct
#' @export
qread.gct <- function(file, type, ...) {
	sep <- "\t";

	if (is.character(file)) {
		f <- base::file(file, "rt");
	} else {
		f <- file;
	}

	if (scan(f, character(), nlines=1, quiet=TRUE) != "#1.2") {
		close(f);
		stop("Input file is not in GCT v1.2 format");
	}

	dims <- scan(f, integer(), nlines=1, sep=sep, quiet=TRUE);
	# dims[1] == number of pobes
	# dims[2] == number of samples

	d <- read.table(f, header=TRUE, sep=sep, quote="", comment.char="", na.strings="",
		colClasses=c("character", "character", rep("numeric", dims[2])),
		nrows = dims[1], stringsAsFactors=FALSE);

	if (is.character(file)) {
		close(f);
	}

	meta <- d[, 1:2];
	data <- d[, 3:ncol(d)];
	rownames(data) <- meta[,1];

	structure(list(meta=meta, data=data), class="annotated.matrix")
}

#' @method qwrite gct
#' @export
qwrite.gct <- function(x, file, type, ...) {
	sep <- "\t";

	if (inherits(x, "annotated.matrix")) {
		data <- x$data;
		meta <- x$meta;
	} else if (is.matrix(x)) {
		data <- x;
		meta <- data.frame(name=rownames(x), description="na");
	} else {
		stop("x is not a matrix or an annotated matrix");
	}

	if (is.character(file)) {
		f <- base::file(file, "wt");
	} else {
		f <- file;
	}

	# write preamble metadata
	writeLines(c("#1.2", paste(as.character(dim(data)), collapse=sep)), f);

	d <- cbind(meta, data);

	write.table(d, f, row.names=FALSE, col.names=TRUE, sep=sep,
		quote=FALSE, na="");

	if (is.character(file)) close(f);
}
