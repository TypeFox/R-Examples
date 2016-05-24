.read_compressed <- function(file, type, ...) {
	if (is.character(file)) {
		file <- as.filename(file);
	}
	if (!is.filename(file)) {
		stop("file must be a filename or a character vector");
	}

	# open conneciton to compressed file, dispatch, and clean up
	f <- .infer_open_method(type)(as.character(file), open="rt");
	x <- qread(f, .infer_file_type(trim_ext(file)), ...);
	close(f);

	x
}

.write_compressed <- function(x, file, type, ...) {
	if (is.character(file)) {
		file <- as.filename(file);
	}
	if (!is.filename(file)) {
		stop("file must be a filename or a character vector");
	}

	# open conneciton to compressed file, dispatch, and clean up
	f <- .infer_open_method(type)(as.character(file), "wt");
	qwrite(x, f, .infer_file_type(trim_ext(file)), ...);
	close(f);
}

.infer_open_method <- function(type) {
	if (type == "gz") {
		gzfile
	} else if (type == "bz") {
		bzfile
	} else if (type == "xz") {
		xzfile
	}
}

#' @method qread gz
#' @export
qread.gz <- .read_compressed;

#' @method qwrite gz
#' @export
qwrite.gz <- .write_compressed;

#' @method qread bz
#' @export
qread.bz <- .read_compressed;

#' @method qwrite bz
#' @export
qwrite.bz <- .write_compressed;

#' @method qread xz
#' @export
qread.xz <- .read_compressed;

#' @method qwrite xz
#' @export
qwrite.xz <- .write_compressed;
