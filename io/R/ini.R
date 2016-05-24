#' @method qread ini
#' @import stringr
#' @export
qread.ini <- function(file, type, ...) {

	inlines <- readLines(file);

	x <- list();
	section <- list();
	section.wrap <- list();

	for (line in inlines) {

		# strip whitespace from ends
		ln <- stringr::str_trim(line);

		if (stringr::str_length(ln) == 0) next;

		# attempt to match section
		m <- stringr::str_match(ln, "^\\[(.+)\\]$")[,2];
		if (!is.na(m)) {
			# finalize and add last section
			if (length(section.wrap) > 0) {
				x <- c(x, section.wrap);
			}

			# create new section
			section.wrap <- list(section);
			names(section.wrap) <- m;

			next
		}
		
		# attempt to match key-value entry
		m <- stringr::str_match(ln, "^(.+)\\s*=\\s*\"?(.+)\"?$");
		if (!is.na(m[,1])) {
			# create key-value pair
			kvp <- list(m[,3]);
			names(kvp) <- list(m[,2]);
			section.wrap[[1]] <- c(section.wrap[[1]], kvp);

			next
		}

		stop("line cannot be parsed: ", ln);

	}

	# finalize and add last section
	if (length(section.wrap) > 0) {
		x <- c(x, section.wrap);
	}

	x
}

#' @method qwrite ini
#' @export
qwrite.ini <- function(x, file, type, ...) {
	
	depth <- .depth_named_list_of_named_list(x);

	if (depth != 1 && depth != 2) {
		stop("`x` must be a named list or a named list of list");
	}

	if (is.character(file)) {
		f <- base::file(file, "wt");
	} else {
		f <- file;
	}

	if (depth == 1) {
		.write_ini_section(f, x)
	} else {
		mapply(
			function(section, x) {
				.write_ini_section(f, x, section)
			},
			names(x),
			x
		);
	}

	if (is.character(file)) {
		close(f);
	}
}


# file is an opened file connection
.write_ini_section <- function(file, named.list, section=NULL) {
	if (!is.null(section)) {
		cat(sprintf("[%s]\n", section), file=file)
	}
	mapply(
		function(key, values) {
			value <- paste(values, collapse=",");
			cat(sprintf("%s=%s\n", key, value), file=file);
		},
		names(named.list),
		named.list
	);
	cat("\n", file=file);
	invisible()
}

.kvdf_to_named_list <- function(d) {
	x <- d[,2];
	names(x) <- d[,1];

	x
}

.depth_named_list_of_named_list <- function(x) {
	depth <- 0;
	z <- x;
	while (is.list(z) && !is.null(names(z))) {
		depth <- depth + 1;
		z <- z[[1]];
	}
	depth
}
