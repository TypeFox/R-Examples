#' Data output
#'
#' This function writes an object to file in a specified format.
#'
#' If \code{type} is \code{NULL}, the file type is inferred from 
#' the file extension. If \code{type} is \code{NA} or if the file extension is
#' unavailable or unknown, \code{type} is inferred from \code{class(x)}.
#' Use \code{\link{io_supported}} to check support for a file or data type.
#' 
#' @param x     data object to write
#' @param file  filename (character or \code{filenamer::filename}),
#'              a readable text-mode connection (for some types),
#'              or path to existing directory
#' @param type  data or file type
#' @param mkpath   whether to create parent directories
#'                 (if they do not already exists)
#' @param symlink  whether to create a symlink to file with a simplified
#'                 file name (ignored if file is not a \code{filename} object);
#'                 an existing file will not be overwritten but an existing
#'                 symlink will be
#' @param ...   other arguments passed to the underlying function
#' @return a data object (object type depends on the underlying function)
#' @export
#' @import filenamer
#'
#' @examples
#' \dontrun{
#' data(cars)
#'
#' # write data to a TSV file
#' qwrite(cars, "cars.tsv")
#' # infer output type based on the class of the cars object
#' qwrite(as.matrix(cars), "cars.mtx", type=NA)
#' }
#'
qwrite <- function(x, file, type=NULL, mkpath=TRUE, symlink=TRUE, ...) {
	if (mkpath && (is.filename(file) || is.character(file))) {
		filenamer::make_path(file);
	}
	.qwrite(
		x,
		if (is.filename(file)) as.character(file) else file,
		if (is.null(type)) .infer_file_type(file) else type,
		...
	);
	if (symlink && (is.filename(file) || is.character(file))) {
		.symlink_simplified(file);
	}

	invisible()
}

.qwrite <- function(x, file, type, ...) {
	# Create a stub variable with the appropriate type
	# prioritize the class of the object to be written
	z <- NA;
	if (is.null(type) || is.na(type)) {
		# infer file type from the class of the object to be written
		class(z) <- class(x);
	} else {
		# use specified/inferred type, with `class(x)` as fallback
		# this allows one to write object in its native text format while
		# using a generic extension (e.g. txt)
		class(z) <- c(type, class(x));
	}

	UseMethod("qwrite", z)
}

#' @method qwrite default
#' @export
qwrite.default <- function(x, file, type, ...) {
	.qio_error(file, type)
}
