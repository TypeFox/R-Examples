#' Data input
#'
#' This function reads a file in a specified format.
#'
#' If \code{type} is \code{NULL}, the file type is inferred from 
#' the file extension.
#' Use \code{\link{io_supported}} to check support for a file or data type.
#' 
#' @param file  file name (character or \code{filenamer::filename}),
#'              a readable text-mode connection (for some types),
#'              or path to existing directory
#' @param type  data or file type
#' @param ...   other arguments passed to the underlying function
#' @return a data object (type depends on the underlying function)
#' @export
#' @import filenamer
#' @import stringr
#'
#' @examples
#' \dontrun{
#' data(cars)
#'
#' # write data to an RDS file
#' qwrite(cars, "cars.rds")
#' # infer output type based on the class of the cars object
#' qwrite(cars, "cars.dfm", type=NA)
#'
#' # read data back in
#' x1 <- qread("cars.rds")
#' # specify the type explicitly
#' x3 <- qread("cars.dfm", type="data.frame")
#' 
#' # read all files (with extension) in current directory
#' xs <- qread(".", pattern="cars")
#' }
#'
qread <- function(file, type=NULL, ...) {
	# Preprocess file name
	if (is.filename(file)) file <- as.character(file);
	
	if (is.character(file)) {
		# Environmental variable substitution
		# NB  The braces {} in the pattern need not be matched... (fix?)
		env.pattern <- "\\$\\{?([a-zA-Z0-9_]+)\\}?";
		locss <- str_locate_all(file, env.pattern);
		var.namess <- lapply(str_match_all(file, env.pattern),
			function(x) {
				if (length(x) > 0) {
					x[, 2]
				} else {
					NULL
				}
			}
		);
		file <- mapply(
			function(f, locs, var.names) {
				# successively substitute environmental variable name with value
				if (nrow(locs) > 0) {
					for (i in 1:nrow(locs)) {
							str_sub(f, locs[i, , drop=FALSE]) <- Sys.getenv(var.names[i])
					}
				}
				f
			},
			file,
			locss,
			var.namess
		);
	}
	
	.qread(
		file,
		if (is.null(type)) .infer_file_type(file) else type,
		...
	)
}

.qread <- function(file, type, ...) {
	# Create a stub variable with the appropriate type
	z <- switch(type,
		matrix = matrix(),
		list = list(),
		# default
		{
			z <- NA; class(z) <- type;
			z
		}
	);
	
	UseMethod("qread", z)
}

#' @method qread default
#' @export
qread.default <- function(file, type, ...) {
	.qio_error(file, type)
}
