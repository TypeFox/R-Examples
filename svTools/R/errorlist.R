### Removes errors concerning a file or a set of files from the list of errors
### Romain Francois <francoisromain@free.fr>
resetErrors <- function (file = NULL, pattern = NULL)
{
	if (!is.null(file)) {
		svTools.env$.errors <- svTools.env$.errors[!svTools.env$.errors$file
			%in% sapply(file, tools::file_path_as_absolute), , drop = FALSE]
	} else if (!is.null(pattern) && any(regexpr(pattern,
		basename(svTools.env$.errors$file)) > 0)) {
		svTools.env$.errors <- svTools.env$.errors[regexpr(pattern,
			svTools.env$.errors$file) <= 0, , drop = FALSE]
	} else { 
		svTools.env$.errors <- emptyError()
	}
	return(invisible(NULL))
}

### Add an error to the list of errors
### Romain Francois <francoisromain@free.fr>
addError <- function (file, line = 1, column = 1, message = "", type = "error")
{
	if (inherits(file, "data.frame") && all(c("file", "line", "column",
		"message", "type") %in% colnames(file))) {
		d <- file
	} else {
		d <- data.frame(file = tools::file_path_as_absolute(file), line = line, 
			column = column, message = message, type = type,
			stringsAsFactors = FALSE)
	}
	svTools.env$.errors <- rbind(svTools.env$.errors, d)
}

### Facility to retrieve the errors stored in the error list
### Romain Francois \<rancoisromain@free.fr>
getErrors <- function (file = NULL, pattern = NULL)
{
	out <- svTools.env$.errors
	if (!is.null(pattern) && any(regexpr(pattern, out$file) > 0))
		out <- out[regexpr(pattern, basename(out$file)) > 0, ]
	if (!is.null(file))
		out <- out[out$file == tools::file_path_as_absolute(file), ]
	if (!is.null(out) && nrow(out)) {
		out <- out[order(out$file, out$line), ]
		out$line <- as.integer(out$line)
		return(out)
	} else return(emptyError()) 
}

### Creates an empty structured error, a data frame with all the columns needed for the
### Romain Francois <francoisromain@free.fr>	
emptyError <- function ()
	return(data.frame(file = character(0), line = integer(0), column = integer(0),
		message = character(0), type = character(0), stringsAsFactors = FALSE))
