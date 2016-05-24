
deparse_to_string <- function (x) {

	if (missing(x)) {
		""
	} else {
		paste0(deparse(x), collapse = "\n")
	}	
}

stopf <- function (string, ...) {
	string <- paste0(string, collapse = '')
	values <- list(...)

	stringified <- lapply(values, function (value) {
		if (!is.character(value)) {
			deparse_to_string(value)
		} else {
			value
		}
	})

	message <- do.call(gettextf, c(list(string), stringified))
	stop(message, call. = FALSE)
}

messagef <- function (string, ...) {
	string <- paste0(string, collapse = '')
	values <- list(...)

	stringified <- lapply(values, function (value) {
		if (!is.character(value)) {
			deparse_to_string(value)
		} else {
			value
		}
	})

	message <- do.call(gettextf, c(list(string), stringified))
	message (message, call. = FALSE)
}

warningf <- function (string, ...) {

	string <- paste0(string, collapse = '')
	values <- list(...)

	stringified <- lapply(values, function (value) {
		if (!is.character(value)) {
			deparse_to_string(value)
		} else {
			value
		}
	})

	message <- do.call(gettextf, c(list(string), stringified))
	warning (message, call. = FALSE)
}

is_boolean <- function (x) {
	is.logical(x) && !is.na(x)
}

