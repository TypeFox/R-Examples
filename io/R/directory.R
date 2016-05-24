#' List the files in a directory.
#'
#' This function extends \code{\link{list.files}} 
#' by excluding the listing of directories.
#' 
#' @param path        a character vector of path names
#' @param full.names  whether to return absolute paths 
#' @param ...         other arguments passed to \code{\link{list.files}}
#' @return a \code{character} vector of only names of files 
#' @export
#'
#' @examples
#' list.files(R.home())
#' list_files(R.home())
#'
list_files <- function(path=".", full.names=FALSE,  ...) {
	x <- list.files(path=path,  full.names=full.names, ...);
	if (!full.names) {
		full.paths <- file.path(path, x);
	} else {
		full.paths <- x;
	}

	x[!file.info(full.paths)$isdir]
}

# Read files from directory
#' @method qread directory
#' @export
qread.directory <- function(file, type, pattern=".*\\..*", closures=FALSE, ...) {
	if (!is.character(file) || !file.info(file)$isdir) {
		stop("`file` should point to a directory with the input files")
	}

	path <- file;

	# Construct file names
	fnames <- list_files(path, full.names=FALSE, pattern=pattern);
	# Prepend path
	fnames.full <- file.path(path, fnames);
	names(fnames.full) <- fnames;

	if (closures) {
		# return list of closures
		lapply(fnames.full,
			function(fn) {
				x <- fn;
				function(...) qread(x, ...)
			}
		)
	} else {
		# return list of data
		lapply(fnames.full, function(fn) qread(fn, ...))
	}
}

# Write objects to directory
#' @method qwrite directory
#' @export
qwrite.directory <- function(x, file, type, file.types=NULL, ...) {
	if (!is.list(x) || is.null(names(x))) {
		stop("`x` must be a named list")
	}
	if (!is.character(file) || !file.info(file)$isdir) {
		stop("`file` should point to a directory where the outputs are to be written")
	}

	if (is.null(file.types)) {
		mapply(
			function(obj, path) {
				qwrite(obj, path, type=NULL, ...)
			},
			x,
			file.path(file, names(x))
		)
	} else {
		mapply(
			function(obj, path, type) {
				qwrite(obj, path, type, ...)
			},
			x,
			file.path(file, names(x)),
			file.types
		)
	}

	invisible()
}

