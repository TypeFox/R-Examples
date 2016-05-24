print.tibbrConnection <-
function(x, ...) {
	cat("Connection to tibbr at ", x$server, " as user ", x$user$login, ".\n", sep="")
	invisible(x)
}

print.tibbrSubject <-
function(x, ...) {
	cat("tibbr subject \"", x$display_name, "\" (", x$name, ").\n", sep="")
	invisible(x)
}

print.tibbrUser <-
function(x, ...) {
	fullname <- if(is.null(x$display_name)) paste(x$first_name, x$last_name) else x$display_name
	cat("tibbr user \"", fullname, "\" (", x$login, ").\n", sep="")
	invisible(x)
}

postingName <- function(x) UseMethod("postingName")

postingName.tibbrSubject <-
function(x) {
	paste("@", as.character(x$name), sep="")
}

postingName.tibbrUser <-
function(x) {
	paste("@", as.character(x$login), sep="")
}

postingName.character <-
function(x) {
	if(regexec("^@?[abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ01234567890_.]+$", x)[[1]] == -1)
		stop("illegal subject name \"", x, "\"", sep="")
	if(substr(x, 0, 1) == "@")
		x
	else
		paste("@", x, sep="")
}

