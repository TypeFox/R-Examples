.infer_file_type <- function(filename) {
	type <- NULL;
	if (is.character(filename)) {
		filename <- as.filename(filename);
	}
	if (is.filename(filename)) {
		if (!is.null(filename$ext)) {
			type <- filename$ext[length(filename$ext)];
		} else {
			# filename has no period: assume it points to a directory
			type <- "directory";
		}
	}
	# `type` may remain NULL if extension is not available
	# or filename actually points to an open connection

	if (is.null(type)) NULL else tolower(type)
}

.qio_error <- function(file, type) {
	if (missing(type)) {
		stop("File type is not specified and could not be inferred from file: ",
			as.character(file));
	} else {
		stop(sprintf("File type `%s` is not supported", type));
	}
}

.check_package <- function(pkg) {
	if (requireNamespace(pkg, quietly=TRUE)) {
		library(pkg, character.only=TRUE, quietly=TRUE);
	} else {
		stop("Please install ", pkg, " or use another file format", call. = FALSE);
	}
}

# To support R < 3.2, provide trimws definition
trimws <- function (x, which = c("both", "left", "right")) {
	which <- match.arg(which)
	mysub <- function(re, x) sub(re, "", x, perl = TRUE)
	if (which == "left") 
		return(mysub("^[ \t\r\n]+", x))
	if (which == "right") 
		return(mysub("[ \t\r\n]+$", x))
	mysub("[ \t\r\n]+$", mysub("^[ \t\r\n]+", x))
}
